version 1.0

import "../../tasks/bioinfo_utils.wdl" as utils

workflow PrepareReads {
    meta {
        description: "## Prepare reads workflow \n Get a sample's reads as FASTQ, whatever they arrived as. Reads that are already FASTQ are passed through, and reads in a CRAM or a BAM are converted. This is the input handling that the mapping workflows all need to do before they can chunk reads and map them; it is a workflow of its own so that a caller that maps the same reads more than once can convert them once and hand the same FASTQs to every run."
    }

    parameter_meta {
        INPUT_READ_FILE_1: "(OPTIONAL) Input sample 1st read pair fastq.gz or fastq"
        INPUT_READ_FILE_2: "(OPTIONAL) Input sample 2nd read pair fastq.gz or fastq"
        INPUT_CRAM_FILE: "(OPTIONAL) Input CRAM file to convert. Needs CRAM_REF and CRAM_REF_INDEX."
        CRAM_REF: "(OPTIONAL) Genome fasta file associated with the CRAM file"
        CRAM_REF_INDEX: "(OPTIONAL) Index of the fasta file associated with the CRAM file"
        INPUT_BAM_FILE: "(OPTIONAL) Input BAM file to convert"
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        INTERLEAVED_READS: "Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'."
        SPLIT_READ_CORES: "Number of cores to use when converting the reads. Default is 8."
        SPLIT_READ_MEM: "Memory, in GB, to use when converting the reads. Default is 50."
    }

    input {
        File? INPUT_READ_FILE_1
        File? INPUT_READ_FILE_2
        File? INPUT_CRAM_FILE
        File? CRAM_REF
        File? CRAM_REF_INDEX
        File? INPUT_BAM_FILE
        Boolean PAIRED_READS = true
        Boolean INTERLEAVED_READS = false
        Int SPLIT_READ_CORES = 8
        Int SPLIT_READ_MEM = 50
    }

    if (defined(INPUT_CRAM_FILE) && defined(CRAM_REF) && defined(CRAM_REF_INDEX)) {
        call utils.convertCRAMtoFASTQ {
            input:
            in_cram_file=INPUT_CRAM_FILE,
            in_ref_file=CRAM_REF,
            in_ref_index_file=CRAM_REF_INDEX,
            in_paired_reads=PAIRED_READS,
            in_cores=SPLIT_READ_CORES,
            in_memory=SPLIT_READ_MEM
        }
    }

    if (defined(INPUT_BAM_FILE)) {
        call utils.convertBAMtoFASTQ {
            input:
            in_bam_file=INPUT_BAM_FILE,
            in_paired_reads=PAIRED_READS,
            in_cores=SPLIT_READ_CORES,
            in_memory=SPLIT_READ_MEM
        }
    }

    File first_read_file = select_first([INPUT_READ_FILE_1, convertCRAMtoFASTQ.output_fastq_1_file, convertBAMtoFASTQ.output_fastq_1_file])

    if (PAIRED_READS && !INTERLEAVED_READS) {
        # There is a real mate file to find, either given to us or made by the
        # conversion.
        File mate_read_file = select_first([INPUT_READ_FILE_2, convertCRAMtoFASTQ.output_fastq_2_file, convertBAMtoFASTQ.output_fastq_2_file])
    }

    output {
        File read_1_file = first_read_file
        # Interleaved reads are all in the first file, so normally there is no
        # second file in that case. Anything the caller passed in as a second
        # FASTQ anyway is still handed back, rather than silently dropped.
        File? read_2_file = if PAIRED_READS && !INTERLEAVED_READS then mate_read_file else INPUT_READ_FILE_2
    }
}
