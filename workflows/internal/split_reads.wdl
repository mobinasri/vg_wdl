version 1.0

import "../../tasks/bioinfo_utils.wdl" as utils

workflow SplitReads {
    meta {
        description: "## Split reads workflow \n Turn a sample's reads into evenly sized FASTQ chunks ready to be mapped in parallel, converting a CRAM or a BAM first if that is what the reads arrived in. The whole, un-split FASTQs can be handed back too, for steps that need to read all of a sample's reads at once, such as counting kmers."
    }

    parameter_meta {
        INPUT_READ_FILE_1: "(OPTIONAL) Input sample 1st read pair fastq.gz or fastq"
        INPUT_READ_FILE_2: "(OPTIONAL) Input sample 2nd read pair fastq.gz or fastq"
        INPUT_CRAM_FILE: "(OPTIONAL) Input CRAM file to convert. Needs CRAM_REF and CRAM_REF_INDEX."
        CRAM_REF: "(OPTIONAL) Genome fasta file associated with the CRAM file"
        CRAM_REF_INDEX: "(OPTIONAL) Index of the fasta file associated with the CRAM file"
        INPUT_BAM_FILE: "(OPTIONAL) Input BAM file to convert"
        OUTPUT_WHOLE_READS: "Should the whole, un-split FASTQs be handed back as well as the chunks? Haplotype sampling needs to see all of a sample's reads at once, so a caller that is going to sample wants these. Default is 'false'."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        INTERLEAVED_READS: "Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'."
        READS_PER_CHUNK: "Number of reads to put in each chunk. Default 20 million."
        SPLIT_READ_CORES: "Number of cores to use when converting and splitting the reads. Default is 8."
        SPLIT_READ_MEM: "Memory, in GB, to use when converting the reads. Default is 50."
    }

    input {
        File? INPUT_READ_FILE_1
        File? INPUT_READ_FILE_2
        File? INPUT_CRAM_FILE
        File? CRAM_REF
        File? CRAM_REF_INDEX
        File? INPUT_BAM_FILE
        Boolean OUTPUT_WHOLE_READS = false
        Boolean PAIRED_READS = true
        Boolean INTERLEAVED_READS = false
        Int READS_PER_CHUNK = 20000000
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
        File second_read_file = select_first([INPUT_READ_FILE_2, convertCRAMtoFASTQ.output_fastq_2_file, convertBAMtoFASTQ.output_fastq_2_file])
    }

    call utils.splitReads as firstReadPair {
        input:
            in_read_file=first_read_file,
            in_pair_id="1",
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES
    }

    if (PAIRED_READS && !INTERLEAVED_READS) {
        call utils.splitReads as secondReadPair {
            input:
            in_read_file=select_first([second_read_file]),
            in_pair_id="2",
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES
        }
    }

    # Only output whole read files if asked
    if (OUTPUT_WHOLE_READS) {
        File whole_first_read_file = first_read_file
        File? whole_second_read_file = second_read_file
    }

    output {
        Array[File] read_chunks_1 = firstReadPair.output_read_chunks
        # Null when the reads are single-ended or interleaved, in which case
        # everything is in the first set of chunks.
        Array[File]? read_chunks_2 = secondReadPair.output_read_chunks
        File? read_1_file = whole_first_read_file
        File? read_2_file = whole_second_read_file
    }
}
