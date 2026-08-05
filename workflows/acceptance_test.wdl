version 1.0

import "../tasks/bioinfo_utils.wdl" as utils
import "./giraffe.wdl" as giraffe_wf
import "./giraffe_and_deepvariant_fromGAF.wdl" as gaf_wf
import "./internal/index_for_giraffe.wdl" as index_wf
import "./internal/prepare_reference.wdl" as reference_wf
import "./internal/split_reads.wdl" as reads_wf

workflow AcceptanceTest {

    meta {
        description: "### Acceptance testing workflow \n Workflow for comparing the calling accuracy for two different vg versions, a \"candidate\" and a \"baseline\". \n\n Runs indexing, mapping, and surjection stages, followed by calling with [DeepVariant](https://github.com/google/deepvariant) and evaluation against a truth set with [Aardvark](https://github.com/PacificBiosciences/aardvark). You can use different vg versions or the same vg version for each stage; stages that don't need to run separately will be run once. By default, indexing is done once, and mapping and surjection are done independently. \n\n Indexes can be provided to skip the indexing stage. If the candidate vg cannot use the baseline vg's indexes, set *CANDIDATE_SEPARATE_INDEXES*. The candidate run then uses the *CANDIDATE_\** index inputs, and anything not passed there is built with the candidate container. With haplotype sampling on, the sampled graph and its indexes count as indexes: they are made once with the baseline vg by default, and once per run when *CANDIDATE_SEPARATE_INDEXES* is set. \n\n More information at [https://github.com/vgteam/vg_wdl/tree/master#giraffe-acceptance-test-workflow](https://github.com/vgteam/vg_wdl/tree/master#acceptance-testing-workflow)."
    }

    parameter_meta {
        BASELINE_VG_DOCKER: "Container image to use when running vg for the baseline run, which is the known-good version to compare against"
        CANDIDATE_VG_DOCKER: "Container image to use when running vg for the candidate run, which is the version under test"
        BASELINE_VG_GIRAFFE_DOCKER: "(OPTIONAL) Container image to use when running vg giraffe mapping in the baseline run, instead of BASELINE_VG_DOCKER"
        BASELINE_VG_SURJECT_DOCKER: "(OPTIONAL) Container image to use when running vg surject in the baseline run, instead of BASELINE_VG_DOCKER"
        CANDIDATE_VG_GIRAFFE_DOCKER: "(OPTIONAL) Container image to use when running vg giraffe mapping in the candidate run, instead of CANDIDATE_VG_DOCKER. Set this and BASELINE_VG_GIRAFFE_DOCKER to the same image to map once and compare only what happens after mapping."
        CANDIDATE_VG_SURJECT_DOCKER: "(OPTIONAL) Container image to use when running vg surject in the candidate run, instead of CANDIDATE_VG_DOCKER"
        BASELINE_VG_SURJECT_OPTIONS: "(OPTIONAL) Extra command line options for vg surject in the baseline run. Set this and CANDIDATE_VG_SURJECT_OPTIONS differently to compare surjection settings rather than, or as well as, vg versions."
        CANDIDATE_VG_SURJECT_OPTIONS: "(OPTIONAL) Extra command line options for vg surject in the candidate run"
        INPUT_READ_FILE_1: "Input sample 1st read pair fastq.gz"
        INPUT_READ_FILE_2: "Input sample 2nd read pair fastq.gz"
        INPUT_CRAM_FILE: "Input CRAM file. Converted to FASTQ once and shared by both runs."
        CRAM_REF: "Genome fasta file associated with the CRAM file"
        CRAM_REF_INDEX: "Index of the fasta file associated with the CRAM file"
        GBZ_FILE: "Path to .gbz index file. Used by both runs unless CANDIDATE_SEPARATE_INDEXES is set."
        DIST_FILE: "(OPTIONAL) Path to .dist index file. Built with the baseline vg if not provided."
        MIN_FILE: "(OPTIONAL) Path to .min index file. Built with the baseline vg if not provided."
        ZIPCODES_FILE: "(OPTIONAL) For chaining-based alignment, path to .zipcodes index file matching MIN_FILE"
        HAPL_FILE: "(OPTIONAL) Path to .hapl file used in haplotype sampling"
        CANDIDATE_SEPARATE_INDEXES: "Should the candidate run get its own indexes instead of sharing the baseline run's? Set this when the two vg versions cannot use each other's indexes. Default is 'false'."
        CANDIDATE_GBZ_FILE: "(OPTIONAL) Path to .gbz index file for the candidate run. Only used if CANDIDATE_SEPARATE_INDEXES is set; defaults to GBZ_FILE."
        CANDIDATE_DIST_FILE: "(OPTIONAL) Path to .dist index file for the candidate run. Only used if CANDIDATE_SEPARATE_INDEXES is set; built with the candidate vg if not provided."
        CANDIDATE_MIN_FILE: "(OPTIONAL) Path to .min index file for the candidate run. Only used if CANDIDATE_SEPARATE_INDEXES is set; built with the candidate vg if not provided."
        CANDIDATE_ZIPCODES_FILE: "(OPTIONAL) Path to .zipcodes index file for the candidate run, matching CANDIDATE_MIN_FILE. Only used if CANDIDATE_SEPARATE_INDEXES is set."
        CANDIDATE_HAPL_FILE: "(OPTIONAL) Path to .hapl file for the candidate run. Only used if CANDIDATE_SEPARATE_INDEXES is set."
        SAMPLE_NAME: "The sample name"
        TRUTH_VCF: "Path to .vcf.gz of truth calls to evaluate both runs against"
        TRUTH_VCF_INDEX: "(OPTIONAL) Tabix index for TRUTH_VCF. Made if not provided."
        EVALUATION_REGIONS_BED: "BED of regions to evaluate in. Required, because Aardvark needs to be told where the truth set is complete."
        STRATIFICATION_ARCHIVE: "(OPTIONAL) tar.gz of a GIAB-style stratification folder (root TSV plus its referenced BED files) to break the Aardvark results down by"
        OUTPUT_GAF: "Should a GAF file with the aligned reads be saved for each run? When both runs map the same way there is only one set of alignments, so both outputs are the same file. Default is 'false'."
        OUTPUT_BAM: "Should the merged BAM be saved for each run? Default is 'false'."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        INTERLEAVED_READS: "Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'."
        READS_PER_CHUNK: "Number of reads contained in each mapping chunk. Default 20 million."
        CONTIGS: "(OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index."
        PATH_LIST_FILE: "(OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths."
        REFERENCE_PREFIX: "Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)"
        REFERENCE_FILE: "(OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths."
        HAPLOID_CONTIGS: "(OPTIONAL) Names of contigs in the reference (without REFERENCE_PREFIX) that are haploid in this sample (often chrX and chrY). Not compatible with DeepVariant 1.5."
        PAR_REGIONS_BED_FILE: "(OPTIONAL) BED file with pseudo-autosomal regions. Not compatible with DeepVariant 1.5."
        PRUNE_LOW_COMPLEXITY: "Whether or not to remove low-complexity or short in-tail anchors when surjecting and force tail realingment. Default is 'true'."
        LEFTALIGN_BAM: "Whether or not to left-align reads in the BAM. Default is 'true'."
        REALIGN_INDELS: "Whether or not to realign reads near indels. Default is 'true'."
        REALIGNMENT_EXPANSION_BASES: "Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160."
        MIN_MAPQ: "Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong. Default is the DeepVariant default for the model type."
        MAX_FRAGMENT_LENGTH: "Maximum distance at which to mark paired reads properly paired. Default is 3000."
        GIRAFFE_PRESET: "(OPTIONAL) Name of Giraffe mapper parameter preset to use (default, fast, hifi, or r10)"
        GIRAFFE_OPTIONS: "(OPTIONAL) Extra command line options for Giraffe mapper"
        DV_MODEL_TYPE: "Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA."
        DV_MODEL_META: ".meta file for a custom DeepVariant calling model"
        DV_MODEL_INDEX: ".index file for a custom DeepVariant calling model"
        DV_MODEL_DATA: ".data-00000-of-00001 file for a custom DeepVariant calling model"
        DV_MODEL_FILES: "Array of all files in the root directory of the DV model, if not using DV_MODEL_META/DV_MODEL_INDEX/DV_MODEL_DATA format"
        DV_MODEL_VARIABLES_FILES: "Array of files that need to go in a 'variables' subdirectory for a DV model"
        DV_KEEP_LEGACY_AC: "Should DV use the legacy allele counter behavior? If unspecified this is not done, unless set in the model. Might want to be on for short reads."
        DV_NORM_READS: "Should DV normalize reads itself? If unspecified this is not done, unless set in the model."
        OTHER_MAKEEXAMPLES_ARG: "Additional arguments for the make_examples step of DeepVariant"
        DV_USE_GPUS: "Should DeepVariant use GPUs for calling variants? Default is 'true'."
        DV_NO_GPU_DOCKER: "Container image to use when running DeepVariant for steps that don't benefit from GPUs. Must be DeepVariant 1.8+."
        DV_GPU_DOCKER: "Container image to use when running DeepVariant for steps that benefit from GPUs. Must be DeepVariant 1.8+."
        SPLIT_READ_CORES: "Number of cores to use when splitting the reads into chunks. Default is 8."
        SPLIT_READ_MEM: "Memory, in GB, to use when splitting the reads into chunks. Default is 50."
        MAP_CORES: "Number of cores to use when mapping the reads. Default is 16."
        MAP_MEM: "Memory, in GB, to use when mapping the reads. Default is 120."
        HAPLOTYPE_SAMPLING: "Whether or not to use haplotype sampling before running giraffe. The sampled graph and its indexes count as indexes, so they are made once unless CANDIDATE_SEPARATE_INDEXES is set. Default is 'true'."
        INDEX_MINIMIZER_WEIGHTED: "Whether to use weighted minimizer indexing. (Default: true)"
        INDEX_MINIMIZER_MEM: "Memory, in GB, to use when making the minimizer index. (Default: 320 if weighted, 120 otherwise)"
        KMER_COUNTING_MEM: "Memory, in GB, to use when counting kmers. (Default: 64)"
        HAPLOTYPE_INDEXING_MEM: "Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling, and giraffe distance index). (Default: 120)"
        BAM_PREPROCESS_MEM: "Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20."
        REALIGN_MEM: "Memory, in GB, to use for Abra indel realignment. Default is 40 or MAP_MEM, whichever is lower."
        CALL_CORES: "Number of cores to use when calling variants. Default is 8."
        CALL_MEM: "Memory, in GB, to use when calling variants. Default is 50."
        MAKE_EXAMPLES_CORES: "Number of cores to use when making DeepVariant examples. Default is CALL_CORES."
        MAKE_EXAMPLES_MEM: "Memory, in GB, to use when making DeepVariant examples. Default is CALL_MEM."
        AARDVARK_CORES: "Number of cores to use when running Aardvark. Default is 16."
        AARDVARK_MEM: "Memory, in GB, to use when running Aardvark. Default is 30."
    }

    input {
        String BASELINE_VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
        String CANDIDATE_VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
        String? BASELINE_VG_GIRAFFE_DOCKER
        String? BASELINE_VG_SURJECT_DOCKER
        String? CANDIDATE_VG_GIRAFFE_DOCKER
        String? CANDIDATE_VG_SURJECT_DOCKER
        String BASELINE_VG_SURJECT_OPTIONS = ""
        String CANDIDATE_VG_SURJECT_OPTIONS = ""
        File? INPUT_READ_FILE_1
        File? INPUT_READ_FILE_2
        File? INPUT_CRAM_FILE
        File? CRAM_REF
        File? CRAM_REF_INDEX
        File GBZ_FILE
        File? DIST_FILE
        File? MIN_FILE
        File? ZIPCODES_FILE
        File? HAPL_FILE
        Boolean CANDIDATE_SEPARATE_INDEXES = false
        File? CANDIDATE_GBZ_FILE
        File? CANDIDATE_DIST_FILE
        File? CANDIDATE_MIN_FILE
        File? CANDIDATE_ZIPCODES_FILE
        File? CANDIDATE_HAPL_FILE
        String SAMPLE_NAME
        File TRUTH_VCF
        File? TRUTH_VCF_INDEX
        File EVALUATION_REGIONS_BED
        File? STRATIFICATION_ARCHIVE
        Boolean OUTPUT_GAF = false
        Boolean OUTPUT_BAM = false
        Boolean PAIRED_READS = true
        Boolean INTERLEAVED_READS = false
        Int READS_PER_CHUNK = 20000000
        Array[String]+? CONTIGS
        File? PATH_LIST_FILE
        String REFERENCE_PREFIX = ""
        File? REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        File? REFERENCE_DICT_FILE
        Array[String]? HAPLOID_CONTIGS
        File? PAR_REGIONS_BED_FILE
        Boolean PRUNE_LOW_COMPLEXITY = true
        Boolean LEFTALIGN_BAM = true
        Boolean REALIGN_INDELS = true
        Int REALIGNMENT_EXPANSION_BASES = 160
        Int? MIN_MAPQ
        Int MAX_FRAGMENT_LENGTH = 3000
        String GIRAFFE_PRESET = "default"
        String GIRAFFE_OPTIONS = ""
        String DV_MODEL_TYPE = "WGS"
        File? DV_MODEL_META
        File? DV_MODEL_INDEX
        File? DV_MODEL_DATA
        Array[File]? DV_MODEL_FILES
        Array[File]? DV_MODEL_VARIABLES_FILES
        Boolean? DV_KEEP_LEGACY_AC
        Boolean? DV_NORM_READS
        String OTHER_MAKEEXAMPLES_ARG = ""
        Boolean DV_USE_GPUS = true
        String? DV_NO_GPU_DOCKER
        String? DV_GPU_DOCKER
        Int SPLIT_READ_CORES = 8
        Int SPLIT_READ_MEM = 50
        Int MAP_CORES = 16
        Int MAP_MEM = 120
        Boolean HAPLOTYPE_SAMPLING = true
        Boolean INDEX_MINIMIZER_WEIGHTED = true
        Int INDEX_MINIMIZER_MEM = if INDEX_MINIMIZER_WEIGHTED then 320 else 120
        Int KMER_COUNTING_MEM = 64
        Int HAPLOTYPE_INDEXING_MEM = 120
        Int BAM_PREPROCESS_MEM = 20
        Int REALIGN_MEM = if MAP_MEM < 40 then MAP_MEM else 40
        Int CALL_CORES = 8
        Int CALL_MEM = 50
        Int MAKE_EXAMPLES_CORES = CALL_CORES
        Int MAKE_EXAMPLES_MEM = CALL_MEM
        Int AARDVARK_CORES = 16
        Int AARDVARK_MEM = 30
    }

    # Which container each vg stage uses. A stage that gets the same container on
    # both sides is not under test.
    String baseline_giraffe_docker = select_first([BASELINE_VG_GIRAFFE_DOCKER, BASELINE_VG_DOCKER])
    String candidate_giraffe_docker = select_first([CANDIDATE_VG_GIRAFFE_DOCKER, CANDIDATE_VG_DOCKER])
    String baseline_surject_docker = select_first([BASELINE_VG_SURJECT_DOCKER, BASELINE_VG_DOCKER])
    String candidate_surject_docker = select_first([CANDIDATE_VG_SURJECT_DOCKER, CANDIDATE_VG_DOCKER])

    # Separate candidate indexes mean a different graph to map against, so the
    # alignments can't be shared no matter what container maps them.
    Boolean share_mapping = !CANDIDATE_SEPARATE_INDEXES && baseline_giraffe_docker == candidate_giraffe_docker

    ####################################################################
    # Everything the two runs are supposed to have in common is set up #
    # here, once, so that the vg container is the only thing that      #
    # differs between them.                                            #
    ####################################################################

    # Both runs have to see exactly the same reads, split exactly the same way, so
    # the reads are converted and chunked once up front.
    call reads_wf.SplitReads {
        input:
        INPUT_READ_FILE_1=INPUT_READ_FILE_1,
        INPUT_READ_FILE_2=INPUT_READ_FILE_2,
        INPUT_CRAM_FILE=INPUT_CRAM_FILE,
        CRAM_REF=CRAM_REF,
        CRAM_REF_INDEX=CRAM_REF_INDEX,
        # Haplotype sampling counts kmers over all of a sample's reads at once.
        OUTPUT_WHOLE_READS=HAPLOTYPE_SAMPLING,
        PAIRED_READS=PAIRED_READS,
        INTERLEAVED_READS=INTERLEAVED_READS,
        READS_PER_CHUNK=READS_PER_CHUNK,
        SPLIT_READ_CORES=SPLIT_READ_CORES,
        SPLIT_READ_MEM=SPLIT_READ_MEM
    }

    # The contig set and the reference come from the graph as given, before any
    # haplotype sampling, and are shared by both runs and by the evaluation.
    # Sharing them is what makes the two call sets comparable at all: they have
    # to be called on the same contigs against the same reference bases. The
    # baseline vg extracts them, since with shared indexes it is the version
    # that made everything else too.
    call reference_wf.PrepareReference {
        input:
        GBZ_FILE=GBZ_FILE,
        CONTIGS=CONTIGS,
        PATH_LIST_FILE=PATH_LIST_FILE,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        REFERENCE_FILE=REFERENCE_FILE,
        REFERENCE_INDEX_FILE=REFERENCE_INDEX_FILE,
        REFERENCE_DICT_FILE=REFERENCE_DICT_FILE,
        EXTRACT_MEM=MAP_MEM,
        VG_DOCKER=BASELINE_VG_DOCKER
    }

    # Both evaluations need the truth set indexed.
    if (!defined(TRUTH_VCF_INDEX)) {
        call utils.indexVcf as indexTruthVcf {
            input:
            in_vcf=TRUTH_VCF
        }
    }
    File truth_vcf_index = select_first([TRUTH_VCF_INDEX, indexTruthVcf.vcf_index_file])

    ###############################################################
    # Indexes: one set for both runs, or one set per run          #
    ###############################################################

    call index_wf.IndexForGiraffe as baselineIndexes {
        input:
        GBZ_FILE=GBZ_FILE,
        DIST_FILE=DIST_FILE,
        MIN_FILE=MIN_FILE,
        ZIPCODES_FILE=ZIPCODES_FILE,
        HAPL_FILE=HAPL_FILE,
        HAPLOTYPE_SAMPLING=HAPLOTYPE_SAMPLING,
        INPUT_READ_FILE_FIRST=SplitReads.read_1_file,
        INPUT_READ_FILE_SECOND=SplitReads.read_2_file,
        GIRAFFE_PRESET=GIRAFFE_PRESET,
        INDEX_MINIMIZER_WEIGHTED=INDEX_MINIMIZER_WEIGHTED,
        INDEX_MINIMIZER_MEM=INDEX_MINIMIZER_MEM,
        KMER_COUNTING_MEM=KMER_COUNTING_MEM,
        HAPLOTYPE_INDEXING_MEM=HAPLOTYPE_INDEXING_MEM,
        CORES=MAP_CORES,
        VG_DOCKER=BASELINE_VG_DOCKER
    }

    if (CANDIDATE_SEPARATE_INDEXES) {
        # The candidate gets its own indexes, from its own inputs where they are
        # given and from its own vg where they are not.
        call index_wf.IndexForGiraffe as candidateIndexes {
            input:
            GBZ_FILE=select_first([CANDIDATE_GBZ_FILE, GBZ_FILE]),
            DIST_FILE=CANDIDATE_DIST_FILE,
            MIN_FILE=CANDIDATE_MIN_FILE,
            ZIPCODES_FILE=CANDIDATE_ZIPCODES_FILE,
            HAPL_FILE=CANDIDATE_HAPL_FILE,
            HAPLOTYPE_SAMPLING=HAPLOTYPE_SAMPLING,
            INPUT_READ_FILE_FIRST=SplitReads.read_1_file,
            INPUT_READ_FILE_SECOND=SplitReads.read_2_file,
            GIRAFFE_PRESET=GIRAFFE_PRESET,
            INDEX_MINIMIZER_WEIGHTED=INDEX_MINIMIZER_WEIGHTED,
            INDEX_MINIMIZER_MEM=INDEX_MINIMIZER_MEM,
            KMER_COUNTING_MEM=KMER_COUNTING_MEM,
            HAPLOTYPE_INDEXING_MEM=HAPLOTYPE_INDEXING_MEM,
            CORES=MAP_CORES,
            VG_DOCKER=CANDIDATE_VG_DOCKER
        }
    }

    File candidate_gbz_file = select_first([candidateIndexes.gbz_file, baselineIndexes.gbz_file])
    File candidate_dist_file = select_first([candidateIndexes.dist_file, baselineIndexes.dist_file])
    File candidate_min_file = select_first([candidateIndexes.min_file, baselineIndexes.min_file])
    # Zipcodes can legitimately be absent, so we take the candidate's own only
    # when the candidate had its own indexes made at all.
    File? candidate_zipcodes_file = if CANDIDATE_SEPARATE_INDEXES then candidateIndexes.zipcodes_file else baselineIndexes.zipcodes_file

    ###############################################################
    # Map: once if both sides map the same way, otherwise twice    #
    ###############################################################

    call giraffe_wf.Giraffe as baselineMapping {
        input:
        READ_CHUNKS_1=SplitReads.read_chunks_1,
        READ_CHUNKS_2=SplitReads.read_chunks_2,
        GBZ_FILE=baselineIndexes.gbz_file,
        DIST_FILE=baselineIndexes.dist_file,
        MIN_FILE=baselineIndexes.min_file,
        ZIPCODES_FILE=baselineIndexes.zipcodes_file,
        SAMPLE_NAME=SAMPLE_NAME,
        # Surjection is a stage of its own here, so this only maps, and hands
        # back the alignment chunks for surjection to pick up.
        OUTPUT_SINGLE_BAM=false,
        OUTPUT_CALLING_BAMS=false,
        OUTPUT_GAF=OUTPUT_GAF,
        OUTPUT_GAF_CHUNKS=true,
        PAIRED_READS=PAIRED_READS,
        INTERLEAVED_READS=INTERLEAVED_READS,
        PATH_LIST_FILE=PrepareReference.path_list_file,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        REFERENCE_FILE=PrepareReference.reference_file,
        REFERENCE_INDEX_FILE=PrepareReference.reference_index_file,
        REFERENCE_DICT_FILE=PrepareReference.reference_dict_file,
        GIRAFFE_PRESET=GIRAFFE_PRESET,
        GIRAFFE_OPTIONS=GIRAFFE_OPTIONS,
        MAP_CORES=MAP_CORES,
        MAP_MEM=MAP_MEM,
        # The indexes are already made, so mapping never samples again.
        HAPLOTYPE_SAMPLING=false,
        VG_DOCKER=BASELINE_VG_DOCKER,
        VG_GIRAFFE_DOCKER=baseline_giraffe_docker
    }

    if (!share_mapping) {
        call giraffe_wf.Giraffe as candidateMapping {
            input:
            READ_CHUNKS_1=SplitReads.read_chunks_1,
            READ_CHUNKS_2=SplitReads.read_chunks_2,
            GBZ_FILE=candidate_gbz_file,
            DIST_FILE=candidate_dist_file,
            MIN_FILE=candidate_min_file,
            ZIPCODES_FILE=candidate_zipcodes_file,
            SAMPLE_NAME=SAMPLE_NAME,
            OUTPUT_SINGLE_BAM=false,
            OUTPUT_CALLING_BAMS=false,
            OUTPUT_GAF=OUTPUT_GAF,
            OUTPUT_GAF_CHUNKS=true,
            PAIRED_READS=PAIRED_READS,
            INTERLEAVED_READS=INTERLEAVED_READS,
            PATH_LIST_FILE=PrepareReference.path_list_file,
            REFERENCE_PREFIX=REFERENCE_PREFIX,
            REFERENCE_FILE=PrepareReference.reference_file,
            REFERENCE_INDEX_FILE=PrepareReference.reference_index_file,
            REFERENCE_DICT_FILE=PrepareReference.reference_dict_file,
            GIRAFFE_PRESET=GIRAFFE_PRESET,
            GIRAFFE_OPTIONS=GIRAFFE_OPTIONS,
            MAP_CORES=MAP_CORES,
            MAP_MEM=MAP_MEM,
            HAPLOTYPE_SAMPLING=false,
            VG_DOCKER=CANDIDATE_VG_DOCKER,
            VG_GIRAFFE_DOCKER=candidate_giraffe_docker
        }
    }

    # We asked for the chunks, so they are there, but the mapping runs hand them
    # back optionally.
    Array[File] baseline_gaf_chunks = select_first([baselineMapping.output_gaf_chunks])
    Array[File] candidate_gaf_chunks = select_first([candidateMapping.output_gaf_chunks, baselineMapping.output_gaf_chunks])

    ###############################################################
    # Surject: twice, since this is a stage that can be tested     #
    ###############################################################

    # Surjection has to use the graph the reads were mapped to, because the
    # alignments name its nodes. No truth set goes in here, so no hap.py runs:
    # the comparison is Aardvark's, below.
    call gaf_wf.GiraffeDeepVariantFromGAF as baselineCalling {
        input:
        GAF_CHUNKS=baseline_gaf_chunks,
        GBZ_FILE=baselineIndexes.gbz_file,
        SAMPLE_NAME=SAMPLE_NAME,
        OUTPUT_SINGLE_BAM=OUTPUT_BAM,
        OUTPUT_CALLING_BAMS=false,
        PAIRED_READS=PAIRED_READS,
        PATH_LIST_FILE=PrepareReference.path_list_file,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        REFERENCE_FILE=PrepareReference.reference_file,
        REFERENCE_INDEX_FILE=PrepareReference.reference_index_file,
        REFERENCE_DICT_FILE=PrepareReference.reference_dict_file,
        HAPLOID_CONTIGS=HAPLOID_CONTIGS,
        PAR_REGIONS_BED_FILE=PAR_REGIONS_BED_FILE,
        PRUNE_LOW_COMPLEXITY=PRUNE_LOW_COMPLEXITY,
        LEFTALIGN_BAM=LEFTALIGN_BAM,
        REALIGN_INDELS=REALIGN_INDELS,
        REALIGNMENT_EXPANSION_BASES=REALIGNMENT_EXPANSION_BASES,
        MIN_MAPQ=MIN_MAPQ,
        MAX_FRAGMENT_LENGTH=MAX_FRAGMENT_LENGTH,
        TRUTH_VCF=TRUTH_VCF,
        TRUTH_VCF_INDEX=truth_vcf_index,
        EVALUATION_REGIONS_BED=EVALUATION_REGIONS_BED,
        EVALUATE_WITH_AARDVARK=true,
        STRATIFICATION_ARCHIVE=STRATIFICATION_ARCHIVE,
        SURJECT_OPTIONS=BASELINE_VG_SURJECT_OPTIONS,
        DV_MODEL_TYPE=DV_MODEL_TYPE,
        DV_MODEL_META=DV_MODEL_META,
        DV_MODEL_INDEX=DV_MODEL_INDEX,
        DV_MODEL_DATA=DV_MODEL_DATA,
        DV_MODEL_FILES=DV_MODEL_FILES,
        DV_MODEL_VARIABLES_FILES=DV_MODEL_VARIABLES_FILES,
        DV_KEEP_LEGACY_AC=DV_KEEP_LEGACY_AC,
        DV_NORM_READS=DV_NORM_READS,
        OTHER_MAKEEXAMPLES_ARG=OTHER_MAKEEXAMPLES_ARG,
        DV_USE_GPUS=DV_USE_GPUS,
        DV_NO_GPU_DOCKER=DV_NO_GPU_DOCKER,
        DV_GPU_DOCKER=DV_GPU_DOCKER,
        VG_CORES=MAP_CORES,
        VG_MEM=MAP_MEM,
        BAM_PREPROCESS_MEM=BAM_PREPROCESS_MEM,
        REALIGN_MEM=REALIGN_MEM,
        CALL_CORES=CALL_CORES,
        CALL_MEM=CALL_MEM,
        MAKE_EXAMPLES_CORES=MAKE_EXAMPLES_CORES,
        MAKE_EXAMPLES_MEM=MAKE_EXAMPLES_MEM,
        EVAL_CORES=AARDVARK_CORES,
        EVAL_MEM=AARDVARK_MEM,
        VG_DOCKER=BASELINE_VG_DOCKER,
        VG_SURJECT_DOCKER=baseline_surject_docker
    }

    call gaf_wf.GiraffeDeepVariantFromGAF as candidateCalling {
        input:
        GAF_CHUNKS=candidate_gaf_chunks,
        # When mapping is shared this is the baseline's graph, which is the one
        # those alignments were made against.
        GBZ_FILE=candidate_gbz_file,
        SAMPLE_NAME=SAMPLE_NAME,
        OUTPUT_SINGLE_BAM=OUTPUT_BAM,
        OUTPUT_CALLING_BAMS=false,
        PAIRED_READS=PAIRED_READS,
        PATH_LIST_FILE=PrepareReference.path_list_file,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        REFERENCE_FILE=PrepareReference.reference_file,
        REFERENCE_INDEX_FILE=PrepareReference.reference_index_file,
        REFERENCE_DICT_FILE=PrepareReference.reference_dict_file,
        HAPLOID_CONTIGS=HAPLOID_CONTIGS,
        PAR_REGIONS_BED_FILE=PAR_REGIONS_BED_FILE,
        PRUNE_LOW_COMPLEXITY=PRUNE_LOW_COMPLEXITY,
        LEFTALIGN_BAM=LEFTALIGN_BAM,
        REALIGN_INDELS=REALIGN_INDELS,
        REALIGNMENT_EXPANSION_BASES=REALIGNMENT_EXPANSION_BASES,
        MIN_MAPQ=MIN_MAPQ,
        MAX_FRAGMENT_LENGTH=MAX_FRAGMENT_LENGTH,
        TRUTH_VCF=TRUTH_VCF,
        TRUTH_VCF_INDEX=truth_vcf_index,
        EVALUATION_REGIONS_BED=EVALUATION_REGIONS_BED,
        EVALUATE_WITH_AARDVARK=true,
        STRATIFICATION_ARCHIVE=STRATIFICATION_ARCHIVE,
        SURJECT_OPTIONS=CANDIDATE_VG_SURJECT_OPTIONS,
        DV_MODEL_TYPE=DV_MODEL_TYPE,
        DV_MODEL_META=DV_MODEL_META,
        DV_MODEL_INDEX=DV_MODEL_INDEX,
        DV_MODEL_DATA=DV_MODEL_DATA,
        DV_MODEL_FILES=DV_MODEL_FILES,
        DV_MODEL_VARIABLES_FILES=DV_MODEL_VARIABLES_FILES,
        DV_KEEP_LEGACY_AC=DV_KEEP_LEGACY_AC,
        DV_NORM_READS=DV_NORM_READS,
        OTHER_MAKEEXAMPLES_ARG=OTHER_MAKEEXAMPLES_ARG,
        DV_USE_GPUS=DV_USE_GPUS,
        DV_NO_GPU_DOCKER=DV_NO_GPU_DOCKER,
        DV_GPU_DOCKER=DV_GPU_DOCKER,
        VG_CORES=MAP_CORES,
        VG_MEM=MAP_MEM,
        BAM_PREPROCESS_MEM=BAM_PREPROCESS_MEM,
        REALIGN_MEM=REALIGN_MEM,
        CALL_CORES=CALL_CORES,
        CALL_MEM=CALL_MEM,
        MAKE_EXAMPLES_CORES=MAKE_EXAMPLES_CORES,
        MAKE_EXAMPLES_MEM=MAKE_EXAMPLES_MEM,
        EVAL_CORES=AARDVARK_CORES,
        EVAL_MEM=AARDVARK_MEM,
        VG_DOCKER=CANDIDATE_VG_DOCKER,
        VG_SURJECT_DOCKER=candidate_surject_docker
    }

    output {
        # Each side is evaluated inside its own calling run, against the same
        # truth set, so the summaries are already there.
        File? baseline_aardvark_summary = baselineCalling.output_aardvark_summary
        File? candidate_aardvark_summary = candidateCalling.output_aardvark_summary
        Array[File]? baseline_aardvark_all_files = baselineCalling.output_aardvark_all_files
        Array[File]? candidate_aardvark_all_files = candidateCalling.output_aardvark_all_files
        File baseline_vcf = baselineCalling.output_vcf
        File baseline_vcf_index = baselineCalling.output_vcf_index
        File candidate_vcf = candidateCalling.output_vcf
        File candidate_vcf_index = candidateCalling.output_vcf_index
        File? baseline_gaf = baselineMapping.output_gaf
        File? candidate_gaf = if share_mapping then baselineMapping.output_gaf else candidateMapping.output_gaf
        File? baseline_bam = baselineCalling.output_bam
        File? baseline_bam_index = baselineCalling.output_bam_index
        File? candidate_bam = candidateCalling.output_bam
        File? candidate_bam_index = candidateCalling.output_bam_index
    }
}
