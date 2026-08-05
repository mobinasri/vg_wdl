version 1.0

import "../../tasks/bioinfo_utils.wdl" as utils
import "../../tasks/gam_gaf_utils.wdl" as gautils
import "../../tasks/vg_map_hts.wdl" as map

workflow MapReads {
    meta {
        description: "## Map reads workflow \n Map FASTQ reads to a pangenome with `vg giraffe`, in parallel over chunks of reads, and hand back the GAF chunks. This is only the mapping: the alignments are still in graph space, and projecting them onto a linear reference is the Surject workflow's job. The chunks are returned unmerged so that a caller can surject them without paying to merge and re-split, and merged as well if it asks. Kept separate from surjection so that one set of alignments can be surjected more than once, with different vg versions or different settings."
    }

    parameter_meta {
        INPUT_READ_FILE_1: "Reads to map, as fastq.gz or fastq. These have to be FASTQ already; PrepareReads gets there from a CRAM or a BAM."
        INPUT_READ_FILE_2: "(OPTIONAL) Second read pair file, if the reads are paired and not interleaved"
        GBZ_FILE: "Path to .gbz index file"
        DIST_FILE: "Path to .dist index file"
        MIN_FILE: "Path to .min index file"
        ZIPCODES_FILE: "(OPTIONAL) For chaining-based alignment, path to .zipcodes index file"
        SAMPLE_NAME: "The sample name, which goes in the read group and sample fields of the alignments"
        OUTPUT_MERGED_GAF: "Should the chunks also be concatenated into one GAF to hand back? The chunks are returned either way. Default is 'false'."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        INTERLEAVED_READS: "Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'."
        READS_PER_CHUNK: "Number of reads contained in each mapping chunk. Default 20 000 000."
        GIRAFFE_PRESET: "Name of Giraffe mapper parameter preset to use (default, fast, hifi, or r10)"
        GIRAFFE_OPTIONS: "Extra command line options for the Giraffe mapper"
        SPLIT_READ_CORES: "Number of cores to use when splitting the reads into chunks. Default is 8."
        MAP_CORES: "Number of cores to use when mapping the reads. Default is 16."
        MAP_MEM: "Memory, in GB, to use when mapping the reads. Default is 120."
        VG_DOCKER: "Container image to use when running vg giraffe"
    }

    input {
        File INPUT_READ_FILE_1
        File? INPUT_READ_FILE_2
        File GBZ_FILE
        File DIST_FILE
        File MIN_FILE
        File? ZIPCODES_FILE
        String SAMPLE_NAME
        Boolean OUTPUT_MERGED_GAF = false
        Boolean PAIRED_READS = true
        Boolean INTERLEAVED_READS = false
        Int READS_PER_CHUNK = 20000000
        String GIRAFFE_PRESET = "default"
        String GIRAFFE_OPTIONS = ""
        Int SPLIT_READ_CORES = 8
        Int MAP_CORES = 16
        Int MAP_MEM = 120
        String VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
    }

    # Split input reads into chunks for parallelized mapping
    call utils.splitReads as firstReadPair {
        input:
            in_read_file=INPUT_READ_FILE_1,
            in_pair_id="1",
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES
    }

    # We always need to pass a full dict file to surjection later, with lengths,
    # because if we pass just path lists and the paths are not completely
    # contained in the graph (like if we're working on GRCh38 paths in a
    # CHM13-based graph), giraffe won't be able to get the path lengths and will
    # crash.
    # TODO: Somehow this problem is supposed to go away if we pull any GRCh38.
    # prefix off the path names by setting REFERENCE_PREFIX and making sure the
    # prefix isn't in the truth set.
    # See <https://github.com/adamnovak/giraffe-dv-wdl/pull/2#issuecomment-955096920>

    if (PAIRED_READS && !INTERLEAVED_READS) {
        call utils.splitReads as secondReadPair {
            input:
            in_read_file=select_first([INPUT_READ_FILE_2]),
            in_pair_id="2",
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES
        }
        Array[Pair[File,File]] read_pair_chunk_files_list = zip(firstReadPair.output_read_chunks, secondReadPair.output_read_chunks)
        scatter (read_pair_chunk_files in read_pair_chunk_files_list) {
            call map.runVGGIRAFFE as runVGGIRAFFE2file {
                input:
                fastq_file_1=read_pair_chunk_files.left,
                fastq_file_2=read_pair_chunk_files.right,
                in_preset=GIRAFFE_PRESET,
                in_giraffe_options=GIRAFFE_OPTIONS,
                in_gbz_file=GBZ_FILE,
                in_dist_file=DIST_FILE,
                in_zipcodes_file=ZIPCODES_FILE,
                in_min_file=MIN_FILE,
                in_sample_name=SAMPLE_NAME,
                nb_cores=MAP_CORES,
                mem_gb=MAP_MEM,
                vg_docker=VG_DOCKER
            }
        }
    }
    # TODO: invent else
    if (!(PAIRED_READS && !INTERLEAVED_READS)) {
        scatter (read_pair_chunk_file in firstReadPair.output_read_chunks) {
            call map.runVGGIRAFFE as runVGGIRAFFE1file {
                input:
                fastq_file_1=read_pair_chunk_file,
                in_interleaved=INTERLEAVED_READS,
                in_preset=GIRAFFE_PRESET,
                in_giraffe_options=GIRAFFE_OPTIONS,
                in_gbz_file=GBZ_FILE,
                in_dist_file=DIST_FILE,
                in_zipcodes_file=ZIPCODES_FILE,
                in_min_file=MIN_FILE,
                in_sample_name=SAMPLE_NAME,
                nb_cores=MAP_CORES,
                mem_gb=MAP_MEM,
                vg_docker=VG_DOCKER
            }
        }
    }

    Array[File] chunks = select_first([runVGGIRAFFE2file.chunk_gaf_file, runVGGIRAFFE1file.chunk_gaf_file])

    if (OUTPUT_MERGED_GAF) {
        call gautils.mergeGAF {
            input:
            in_sample_name=SAMPLE_NAME,
            in_gaf_chunk_files=chunks
        }
    }

    output {
        Array[File] gaf_chunks = chunks
        File? merged_gaf = mergeGAF.output_merged_gaf
    }
}
