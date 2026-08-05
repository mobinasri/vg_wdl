version 1.0

import "../../tasks/bioinfo_utils.wdl" as utils
import "../../tasks/gam_gaf_utils.wdl" as gautils

workflow Surject {
    meta {
        description: "## Surject workflow \n Project graph alignments onto the reference paths with `vg surject` and hand back one sorted BAM. The GAF chunks are surjected in parallel, then sorted and merged, so this takes the chunks MapReads produces rather than a single merged GAF. No BAM post-processing (left-shifting, indel realignment, splitting by contig) happens here; that belongs to whatever calls variants next. Kept separate from mapping so that one set of alignments can be surjected more than once, with different vg versions or different settings."
    }

    parameter_meta {
        GAF_CHUNKS: "GAF chunks to surject in parallel, as produced by MapReads"
        GBZ_FILE: "Path to the .gbz index file the alignments are against. Has to be the graph the reads were mapped to, since the alignments name its nodes."
        PATH_LIST_FILE: "Text file where each line is a path name in the GBZ index to surject onto"
        REFERENCE_DICT_FILE: "Sequence dictionary of the reference, used to order the BAM's contigs"
        REFERENCE_PREFIX: "Remove this off the beginning of path names in the surjected BAM (set to match prefix in PATH_LIST_FILE)"
        SAMPLE_NAME: "The sample name, which goes in the BAM's read group and sample fields"
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        PRUNE_LOW_COMPLEXITY: "Whether or not to remove low-complexity or short in-tail anchors when surjecting and force tail realingment. Default is 'true'."
        MAX_FRAGMENT_LENGTH: "Maximum distance at which to mark paired reads properly paired. Default is 3000."
        SURJECT_CORES: "Number of cores to use when surjecting. Default is 16."
        SURJECT_MEM: "Memory, in GB, to use when surjecting. Default is 120."
        BAM_PREPROCESS_MEM: "Memory, in GB, to use when sorting and merging BAMs. Default is 20."
        VG_DOCKER: "Container image to use when running vg surject"
    }

    input {
        Array[File] GAF_CHUNKS
        File GBZ_FILE
        File PATH_LIST_FILE
        File REFERENCE_DICT_FILE
        String REFERENCE_PREFIX = ""
        String SAMPLE_NAME
        Boolean PAIRED_READS = true
        Boolean PRUNE_LOW_COMPLEXITY = true
        Int MAX_FRAGMENT_LENGTH = 3000
        Int SURJECT_CORES = 16
        Int SURJECT_MEM = 120
        Int BAM_PREPROCESS_MEM = 20
        String VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
    }

    scatter (gaf_file in GAF_CHUNKS) {
        call gautils.surjectGAFtoBAM {
            input:
            in_gaf_file=gaf_file,
            in_gbz_file=GBZ_FILE,
            in_path_list_file=PATH_LIST_FILE,
            in_sample_name=SAMPLE_NAME,
            in_max_fragment_length=MAX_FRAGMENT_LENGTH,
            in_paired_reads=PAIRED_READS,
            in_prune_low_complexity=PRUNE_LOW_COMPLEXITY,
            nb_cores=SURJECT_CORES,
            mem_gb=SURJECT_MEM,
            vg_docker=VG_DOCKER
        }

        call utils.sortBAM {
            input:
            in_bam_file=surjectGAFtoBAM.output_bam_file,
            in_ref_dict=REFERENCE_DICT_FILE,
            in_prefix_to_strip=REFERENCE_PREFIX,
            nb_cores=SURJECT_CORES,
            mem_gb=BAM_PREPROCESS_MEM
        }
    }

    # Merge up the unprocessed surjected alignments
    call utils.mergeAlignmentBAMChunks {
        input:
        in_sample_name=SAMPLE_NAME,
        in_alignment_bam_chunk_files=sortBAM.sorted_bam,
        in_cores=SURJECT_CORES,
        mem_gb=BAM_PREPROCESS_MEM
    }

    output {
        File bam_file = mergeAlignmentBAMChunks.merged_bam_file
        File bam_index_file = mergeAlignmentBAMChunks.merged_bam_file_index
    }
}
