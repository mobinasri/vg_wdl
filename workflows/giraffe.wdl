version 1.0

import "../tasks/validation.wdl" as validation
import "../tasks/bioinfo_utils.wdl" as utils
import "./haplotype_sampling.wdl" as hapl
import "./internal/map_reads.wdl" as map_wf
import "./internal/prepare_reads.wdl" as reads_wf
import "./internal/prepare_reference.wdl" as reference_wf
import "./internal/surject.wdl" as surject_wf

workflow Giraffe {
    meta {
        description: "## Giraffe workflow \n Core VG Giraffe mapping, usable for DeepVariant. Reads are mapped to a pangenome with vg giraffe and pre-processed (e.g. indel realignment). More information at [https://github.com/vgteam/vg_wdl/tree/master#giraffe-workflow](https://github.com/vgteam/vg_wdl/tree/master#giraffe-workflow)."
    }
    parameter_meta {
        INPUT_READ_FILE_1: "Input sample 1st read pair fastq.gz or fastq"
        INPUT_READ_FILE_2: "Input sample 2nd read pair fastq.gz or fastq"
        INPUT_CRAM_FILE: "Input CRAM file to realign"
        CRAM_REF: "Genome fasta file associated with the CRAM file"
        CRAM_REF_INDEX: "Index of the fasta file associated with the CRAM file"
        INPUT_BAM_FILE: "Input BAM file to realign"
        GBZ_FILE: "Path to .gbz index file"
        DIST_FILE: "Path to .dist index file. Optional if using haplotype sampling."
        MIN_FILE: "Path to .min index file. Optional if using haplotype sampling."
        ZIPCODES_FILE: "(OPTIONAL) For chaining-based alignment, path to .zipcodes index file"
        SAMPLE_NAME: "The sample name"
        OUTPUT_SINGLE_BAM: "Should a single merged BAM file be saved? Default is 'true'."
        OUTPUT_CALLING_BAMS: "Should individual contig BAMs be saved? Default is 'false'."
        OUTPUT_GAF: "Should a GAF file with the aligned reads be saved? Default is 'false'."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        INTERLEAVED_READS: "Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'."
        READS_PER_CHUNK: "Number of reads contained in each mapping chunk. Default 20 000 000."
        PATH_LIST_FILE: "(OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths. If using REFERENCE_PREFIX, contig names in here should have the prefix."
        CONTIGS: "(OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index. If using REFERENCE_PREFIX, contig names in here should have the prefix."
        REFERENCE_PREFIX: "Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)"
        REFERENCE_FILE: "(OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference. If using REFERENCE_PREFIX, contig names in here should not have the prefix."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file. If using REFERENCE_PREFIX, contig names in here should not have the prefix."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if REFERENCE_INDEX_FILE is set. If using REFERENCE_PREFIX, contig names in here should not have the prefix. This is used in BAM processing and not for choosing contigs for the surjection, which uses PATH_LIST_FILE."
        PRUNE_LOW_COMPLEXITY: "Whether or not to remove low-complexity or short in-tail anchors when surjecting and force tail realingment. Default is 'true'."  
        LEFTALIGN_BAM: "Whether or not to left-align reads in the BAM. Default is 'true'."
        REALIGN_INDELS: "Whether or not to realign reads near indels. Default is 'true'."
        REALIGNMENT_EXPANSION_BASES: "Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160."
        MAX_FRAGMENT_LENGTH: "Maximum distance at which to mark paired reads properly paired. Default is 3000."
        GIRAFFE_PRESET: "(OPTIONAL) Name of Giraffe mapper parameter preset to use (default, fast, hifi, or r10)"  
        GIRAFFE_OPTIONS: "(OPTIONAL) extra command line options for Giraffe mapper"
        SPLIT_READ_CORES: "Number of cores to use when splitting the reads into chunks. Default is 8."
        SPLIT_READ_MEM: "Memory, in GB, to use when splitting the reads into chunks. Default is 50."
        MAP_CORES: "Number of cores to use when mapping the reads. Default is 16."
        MAP_MEM: "Memory, in GB, to use when mapping the reads. Default is 120."
        BAM_PREPROCESS_MEM: "Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20."
        REALIGN_MEM: "Memory, in GB, to use for Abra indel realignment. Default is 40 or MAP_MEM, whichever is lower."
        HAPLOTYPE_SAMPLING: "Whether or not to use haplotype sampling before running giraffe. Default is 'true'"
        DIPLOID:"Whether or not to use diploid sampling while doing haplotype sampling. Has to use with Haplotype_sampling=true. Default is 'true'"
        SET_REFERENCE:"(OPTIONAL) Name of the single reference to keep for haplotype sampling."
        HAPL_FILE: "(OPTIONAL) Path to .hapl file used in haplotype sampling"
        R_INDEX_FILE: "(OPTIONAL) Path to .ri file used in haplotype sampling"
        KFF_FILE: "(OPTIONAL) Path to .kff file used in haplotype sampling"
        HAPLOTYPE_NUMBER: "Number of generated synthetic haplotypes used in haplotype sampling. (Default: 32)"
        INDEX_MINIMIZER_WEIGHTED: "Whether to use weighted minimizer indexing with haplotype sampling. (Default: true)"
        INDEX_MINIMIZER_MEM: "Memory, in GB, to use when making the minimizer index. (Default: 320 if weighted, 120 otherwise)"
        KMER_COUNTING_MEM: "Memory, in GB, to use when counting kmers. (Default: 64)"
        HAPLOTYPE_INDEXING_MEM: "Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling, and giraffe distance index). (Default: 120)"

        VG_DOCKER: "Container image to use when running vg"
        VG_GIRAFFE_DOCKER: "Alternate container image to use when running vg giraffe mapping"
        VG_SURJECT_DOCKER: "Alternate container image to use when running vg surject"
    }
    input {
        File? INPUT_READ_FILE_1
        File? INPUT_READ_FILE_2
        File? INPUT_CRAM_FILE
        File? CRAM_REF
        File? CRAM_REF_INDEX
        File? INPUT_BAM_FILE
        File GBZ_FILE
        File? DIST_FILE
        File? MIN_FILE
        File? ZIPCODES_FILE
        String SAMPLE_NAME
        Boolean OUTPUT_SINGLE_BAM = true
        Boolean OUTPUT_CALLING_BAMS = false
        Boolean OUTPUT_GAF = false
        Boolean PAIRED_READS = true
        Boolean INTERLEAVED_READS = false
        Int READS_PER_CHUNK = 20000000
        File? PATH_LIST_FILE
        Array[String]+? CONTIGS
        String REFERENCE_PREFIX = ""
        File? REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        File? REFERENCE_DICT_FILE
        Boolean PRUNE_LOW_COMPLEXITY = true
        Boolean LEFTALIGN_BAM = true
        Boolean REALIGN_INDELS = true
        Int REALIGNMENT_EXPANSION_BASES = 160
        Int MAX_FRAGMENT_LENGTH = 3000
        String GIRAFFE_PRESET = "default"
        String GIRAFFE_OPTIONS = ""
        Int SPLIT_READ_CORES = 8
        Int SPLIT_READ_MEM = 50
        Int MAP_CORES = 16
        Int MAP_MEM = 120
        Int BAM_PREPROCESS_MEM = 20
        Int REALIGN_MEM = if MAP_MEM < 40 then MAP_MEM else 40
        Boolean HAPLOTYPE_SAMPLING = true
        Boolean DIPLOID = true
        String? SET_REFERENCE
        File? HAPL_FILE
        File? R_INDEX_FILE
        File? KFF_FILE
        Int HAPLOTYPE_NUMBER = 32
        Boolean INDEX_MINIMIZER_WEIGHTED = true
        Int INDEX_MINIMIZER_MEM = if INDEX_MINIMIZER_WEIGHTED then 320 else 120
        Int KMER_COUNTING_MEM = 64
        Int HAPLOTYPE_INDEXING_MEM = 120

        String VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
        String? VG_GIRAFFE_DOCKER
        String? VG_SURJECT_DOCKER
    }

    # Set up enough path list stuff to do validation as soon as possible.
    # We don't just check the final path list because that could come from
    # haplotype sampling and we want to check any provided path list before going
    # into haplotype sampling.
    if (defined(CONTIGS)) {
        # Put the paths in a file to use later. We know the value is defined,
        # but WDL is a bit low on unboxing calls for optionals so we use
        # select_first.
        File written_path_names_file = write_lines(select_first([CONTIGS]))
        if (REFERENCE_PREFIX != "") {
            call validation.checkPathList as checkWrittenPathList {
                input:
                    in_path_list_file=written_path_names_file,
                    in_reference_prefix=REFERENCE_PREFIX
            }
        }
    }
    if (defined(PATH_LIST_FILE) && REFERENCE_PREFIX != "") {
        call validation.checkPathList as checkProvidedPathList {
            input:
                in_path_list_file=select_first([PATH_LIST_FILE]),
                in_reference_prefix=REFERENCE_PREFIX
        }
    }

    # Validate the dict file if fed in.
    if (defined(REFERENCE_DICT_FILE) && REFERENCE_PREFIX != "") {
        call validation.checkDict as checkProvidedDict {
            input:
                in_dict_file=select_first([REFERENCE_DICT_FILE]),
                in_reference_prefix=REFERENCE_PREFIX
        }
    }


    # Get the reads as FASTQ, whatever container they came in.
    call reads_wf.PrepareReads {
        input:
        INPUT_READ_FILE_1=INPUT_READ_FILE_1,
        INPUT_READ_FILE_2=INPUT_READ_FILE_2,
        INPUT_CRAM_FILE=INPUT_CRAM_FILE,
        CRAM_REF=CRAM_REF,
        CRAM_REF_INDEX=CRAM_REF_INDEX,
        INPUT_BAM_FILE=INPUT_BAM_FILE,
        PAIRED_READS=PAIRED_READS,
        INTERLEAVED_READS=INTERLEAVED_READS,
        SPLIT_READ_CORES=SPLIT_READ_CORES,
        SPLIT_READ_MEM=SPLIT_READ_MEM
    }
    File read_1_file = PrepareReads.read_1_file
    # Null unless there is a separate mate file, which hap sampling and the
    # paired mapping path both need.
    File? read_2_file = PrepareReads.read_2_file

    if (HAPLOTYPE_SAMPLING) {
        call hapl.HaplotypeSampling {
        input:
            GBZ_FILE=GBZ_FILE,
            INPUT_READ_FILE_FIRST=read_1_file,
            # If we're not doing paired reads the result here is probably null.
            INPUT_READ_FILE_SECOND=read_2_file,
            HAPL_FILE=HAPL_FILE,
            DIST_FILE=DIST_FILE,
            R_INDEX_FILE=R_INDEX_FILE,
            KFF_FILE=KFF_FILE,
            HAPLOTYPE_NUMBER=HAPLOTYPE_NUMBER,
            DIPLOID=DIPLOID,
            SET_REFERENCE=SET_REFERENCE,
            INDEX_MINIMIZER_K = if GIRAFFE_PRESET == "default" || GIRAFFE_PRESET == "fast" then 29 else 31,
            INDEX_MINIMIZER_W = if GIRAFFE_PRESET == "default" || GIRAFFE_PRESET == "fast" then 11 else 50,
            INDEX_MINIMIZER_WEIGHTED=INDEX_MINIMIZER_WEIGHTED,
            CORES=MAP_CORES,
            KMER_COUNTING_MEM=KMER_COUNTING_MEM,
            HAPLOTYPE_INDEXING_MEM=HAPLOTYPE_INDEXING_MEM,
            INDEX_MINIMIZER_MEM=INDEX_MINIMIZER_MEM,
            VG_DOCKER=VG_DOCKER
        }

    }

    File file_gbz = select_first([HaplotypeSampling.sampled_graph, GBZ_FILE])
    File file_min = select_first([HaplotypeSampling.sampled_min, MIN_FILE])
    # The zipcode file is optional but we still have a priority list of places to get it from.
    # But we can't select_first since they all might be null.
    Array[File] possible_zipcode_files = select_all([HaplotypeSampling.sampled_zipcodes, ZIPCODES_FILE])
    # We can't actually use None in WDL 1.0 so we need to use a nonexistent null file.
    if (false) {
        Array[File] no_files = []
        #@ except: UnnecessaryFunctionCall
        File NULL_FILE = select_first(no_files)
    }
    File? file_zipcodes = if length(possible_zipcode_files) > 0 then possible_zipcode_files[0] else NULL_FILE
    File file_dist = select_first([HaplotypeSampling.sampled_dist, DIST_FILE])


    # Which path names to work on, and what reference to surject against? These
    # come from the graph we are actually mapping to, which is the sampled one
    # if we sampled.
    call reference_wf.PrepareReference {
        input:
        GBZ_FILE=file_gbz,
        CONTIGS=CONTIGS,
        PATH_LIST_FILE=PATH_LIST_FILE,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        REFERENCE_FILE=REFERENCE_FILE,
        REFERENCE_INDEX_FILE=REFERENCE_INDEX_FILE,
        REFERENCE_DICT_FILE=REFERENCE_DICT_FILE,
        EXTRACT_MEM=MAP_MEM,
        VG_DOCKER=VG_DOCKER
    }
    File pipeline_path_list_file = PrepareReference.path_list_file
    File reference_file = PrepareReference.reference_file
    File reference_index_file = PrepareReference.reference_index_file
    File reference_dict_file = PrepareReference.reference_dict_file

    ################################################################
    # Distribute vg mapping operation over each chunked read pair #
    ################################################################
    call map_wf.MapReads {
        input:
        INPUT_READ_FILE_1=read_1_file,
        INPUT_READ_FILE_2=read_2_file,
        GBZ_FILE=file_gbz,
        DIST_FILE=file_dist,
        MIN_FILE=file_min,
        ZIPCODES_FILE=file_zipcodes,
        SAMPLE_NAME=SAMPLE_NAME,
        OUTPUT_MERGED_GAF=OUTPUT_GAF,
        PAIRED_READS=PAIRED_READS,
        INTERLEAVED_READS=INTERLEAVED_READS,
        READS_PER_CHUNK=READS_PER_CHUNK,
        GIRAFFE_PRESET=GIRAFFE_PRESET,
        GIRAFFE_OPTIONS=GIRAFFE_OPTIONS,
        SPLIT_READ_CORES=SPLIT_READ_CORES,
        MAP_CORES=MAP_CORES,
        MAP_MEM=MAP_MEM,
        VG_DOCKER=select_first([VG_GIRAFFE_DOCKER, VG_DOCKER])
    }

    if (OUTPUT_SINGLE_BAM || OUTPUT_CALLING_BAMS) {
        # We are outputting BAM so surjection is needed
        call surject_wf.Surject {
            input:
            GAF_CHUNKS=MapReads.gaf_chunks,
            GBZ_FILE=file_gbz,
            PATH_LIST_FILE=pipeline_path_list_file,
            REFERENCE_DICT_FILE=reference_dict_file,
            REFERENCE_PREFIX=REFERENCE_PREFIX,
            SAMPLE_NAME=SAMPLE_NAME,
            PAIRED_READS=PAIRED_READS,
            PRUNE_LOW_COMPLEXITY=PRUNE_LOW_COMPLEXITY,
            MAX_FRAGMENT_LENGTH=MAX_FRAGMENT_LENGTH,
            SURJECT_CORES=MAP_CORES,
            SURJECT_MEM=MAP_MEM,
            BAM_PREPROCESS_MEM=BAM_PREPROCESS_MEM,
            VG_DOCKER=select_first([VG_SURJECT_DOCKER, VG_DOCKER])
        }

        if (OUTPUT_CALLING_BAMS || LEFTALIGN_BAM || REALIGN_INDELS) {
            # We will need to split up the BAM by contig to do processing on it.
            # TODO: Unify BAM processing with deepvariant.wdl

            # Split merged alignment by contigs list
            call utils.splitBAMbyPath {
                input:
                in_sample_name=SAMPLE_NAME,
                in_merged_bam_file=Surject.bam_file,
                in_merged_bam_file_index=Surject.bam_index_file,
                in_path_list_file=pipeline_path_list_file,
                in_prefix_to_strip=REFERENCE_PREFIX,
                thread_count=SPLIT_READ_CORES,
                mem_gb=SPLIT_READ_MEM
            }

            ##
            ## Prepare each BAM
            ##
            scatter (bam_and_index_for_path in zip(splitBAMbyPath.bam_contig_files, splitBAMbyPath.bam_contig_files_index)) {
                ## Evantually shift and realign reads
                if (LEFTALIGN_BAM){
                    # Just left-shift each read individually
                    call utils.leftShiftBAMFile {
                        input:
                        in_bam_file=bam_and_index_for_path.left,
                        in_reference_file=reference_file,
                        in_reference_index_file=reference_index_file,
                        mem_gb=BAM_PREPROCESS_MEM
                    }
                }
                if (REALIGN_INDELS) {
                    File forrealign_bam = select_first([leftShiftBAMFile.output_bam_file, bam_and_index_for_path.left])
                    File forrealign_index = select_first([leftShiftBAMFile.output_bam_index_file, bam_and_index_for_path.right])
                    # Do indel realignment
                    call utils.prepareRealignTargets {
                        input:
                        in_bam_file=forrealign_bam,
                        in_bam_index_file=forrealign_index,
                        in_reference_file=reference_file,
                        in_reference_index_file=reference_index_file,
                        in_reference_dict_file=reference_dict_file,
                        in_expansion_bases=REALIGNMENT_EXPANSION_BASES,
                        thread_count=MAP_CORES,
                        mem_gb=BAM_PREPROCESS_MEM
                    }
                    call utils.runAbraRealigner {
                        input:
                            in_bam_file=forrealign_bam,
                            in_bam_index_file=forrealign_index,
                            in_target_bed_file=prepareRealignTargets.output_target_bed_file,
                            in_reference_file=reference_file,
                            in_reference_index_file=reference_index_file,
                            threadCount=MAP_CORES,
                            memoryGb=REALIGN_MEM
                    }
                }
                File processed_bam = select_first([runAbraRealigner.indel_realigned_bam, leftShiftBAMFile.output_bam_file, bam_and_index_for_path.left])
                File processed_bam_index = select_first([runAbraRealigner.indel_realigned_bam_index, leftShiftBAMFile.output_bam_index_file, bam_and_index_for_path.right])
            }
        
            if (OUTPUT_SINGLE_BAM && (LEFTALIGN_BAM || REALIGN_INDELS)){
                # We're outputting one big BAM and we've actually made changes to the chunks. Put them back together.
                call utils.mergeAlignmentBAMChunks as mergeBAM {
                    input:
                    in_sample_name=SAMPLE_NAME,
                    in_alignment_bam_chunk_files=flatten([processed_bam, [splitBAMbyPath.bam_unmapped_file]]),
                    in_cores=SPLIT_READ_CORES,
                    mem_gb=SPLIT_READ_MEM
                }
            }
            
            if (OUTPUT_CALLING_BAMS){
                Array[File] calling_bams = processed_bam
                Array[File] calling_bam_indexes = processed_bam_index
            }
        }

        if (OUTPUT_SINGLE_BAM) {
            # Find the single BAM and index that we want to output.
            # We want the one after postprocessing if we did any, and the plain merged sorted BAM otherwise.
            File single_bam = select_first([mergeBAM.merged_bam_file, Surject.bam_file])
            File single_bam_index = select_first([mergeBAM.merged_bam_file_index, Surject.bam_index_file])
        }
    }

    output {
        File? output_bam = single_bam
        File? output_bam_index = single_bam_index
        File? output_gaf = MapReads.merged_gaf
        Array[File]? output_calling_bams = calling_bams
        Array[File]? output_calling_bam_indexes = calling_bam_indexes

    }   
}

