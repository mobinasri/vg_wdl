version 1.0

import "../tasks/aardvark_evaluation.wdl" as eval
import "../tasks/bioinfo_utils.wdl" as utils

workflow AardvarkEvaluation {

    meta {
        description: "Evaluate small variants using Aardvark. Follows the same pattern as happy_evaluation.wdl."
    }

    parameter_meta {
        QUERY_VCF: "bgzipped VCF with variant calls to evaluate"
        QUERY_VCF_INDEX: "(Optional) tabix index for QUERY_VCF; will be indexed if not provided"
        TRUTH_VCF: "bgzipped VCF with truthset"
        TRUTH_VCF_INDEX: "(Optional) tabix index for TRUTH_VCF; will be indexed if not provided"
        REFERENCE_FILE: "FASTA reference"
        REFERENCE_INDEX_FILE: "(Optional) .fai index; will be indexed if not provided"
        REGIONS_BED: "BED of regions to restrict comparison to"
        STRATIFICATION_ARCHIVE: "(Optional) tar.gz of a GIAB-style stratification folder (root TSV plus its referenced BED files)"
        SAMPLE_NAME: "Sample name, used to name the output directory/archive"
        THREADS: "Number of threads for aardvark compare. Default 16."
        EVAL_MEM: "Memory, in GB, to use when evaluating variant calls. Default is 30."
    }

    input {
        File QUERY_VCF
        File? QUERY_VCF_INDEX
        File TRUTH_VCF
        File? TRUTH_VCF_INDEX
        File REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        File REGIONS_BED
        File? STRATIFICATION_ARCHIVE
        String SAMPLE_NAME
        Int THREADS = 16
        Int EVAL_MEM = 30
    }



    ## index the reference FASTA if needed
    if (!defined(REFERENCE_INDEX_FILE)) {
        call utils.indexReference {
            input:
            in_reference_file=REFERENCE_FILE
        }
    }
    File reference_index_file = select_first([REFERENCE_INDEX_FILE, indexReference.reference_index_file])

    ## index the query VCF if needed
    if (!defined(QUERY_VCF_INDEX)) {
        call utils.indexVcf {
            input:
            in_vcf=QUERY_VCF
        }
    }
    File query_vcf_index = select_first([QUERY_VCF_INDEX, indexVcf.vcf_index_file])

    ## index the truth VCF if needed
    if (!defined(TRUTH_VCF_INDEX)) {
        call utils.indexVcf as indexTruthVcf {
            input:
            in_vcf=TRUTH_VCF
        }
    }
    File truth_vcf_index = select_first([TRUTH_VCF_INDEX, indexTruthVcf.vcf_index_file])

    call eval.compareCallsAardvark {
        input:
        in_query_vcf_file=QUERY_VCF,
        in_stratification_archive=STRATIFICATION_ARCHIVE,
        in_query_vcf_index_file=query_vcf_index,
        in_truth_vcf_file=TRUTH_VCF,
        in_truth_vcf_index_file=truth_vcf_index,
        in_reference_file=REFERENCE_FILE,
        in_reference_index_file=reference_index_file,
        in_regions_bed=REGIONS_BED,
        in_output_prefix=SAMPLE_NAME,
        in_threads=THREADS,
        in_mem=EVAL_MEM
    }

    output {
        File aardvark_summary = compareCallsAardvark.output_summary
        Array[File] aardvark_all_files = compareCallsAardvark.output_all_files
    }
}