version 1.0


task compareCallsAardvark {
    input {
        File in_query_vcf_file
        File? in_query_vcf_index_file
        File in_truth_vcf_file
        File? in_truth_vcf_index_file
        File in_reference_file
        File? in_reference_index_file
        File in_regions_bed
        File? in_stratification_archive
        String in_output_prefix
        Int in_threads = 16
        Int in_disk = 3 * round(size(in_query_vcf_file, "G") + size(in_truth_vcf_file, "G") + size(in_reference_file, "G")) + 20
        Int in_mem = 30
    }
    command <<<
        set -eux -o pipefail


        # Put files where aardvark expects to find index files alongside them
        ln -s "~{in_query_vcf_file}" query.vcf.gz
        ln -s "~{in_truth_vcf_file}" truth.vcf.gz
        ln -s "~{in_reference_file}" reference.fa
        ~{"ln -s " + in_query_vcf_index_file + " query.vcf.gz.tbi"}
        ~{"ln -s " + in_truth_vcf_index_file + " truth.vcf.gz.tbi"}
        ~{"ln -s " + in_reference_index_file + " reference.fa.fai"}

        STRAT_ARG=""
        if [ -n "~{if defined(in_stratification_archive) then "yes" else ""}" ]; then
            mkdir -p strat_dir
            tar -xzf "~{in_stratification_archive}" -C strat_dir
            STRAT_TSV=$(find strat_dir -name "*.tsv" | head -n1)
            STRAT_ARG="--stratification ${STRAT_TSV}"
        fi

        mkdir -p "~{in_output_prefix}"

        aardvark compare \
            --reference reference.fa \
            --truth-vcf truth.vcf.gz \
            --query-vcf query.vcf.gz \
            --regions "~{in_regions_bed}" \
            --output-dir "~{in_output_prefix}" \
            ${STRAT_ARG} \
            --threads ~{in_threads}


    >>>
    output {
        File output_summary = "~{in_output_prefix}/summary.tsv"
        Array[File] output_all_files = glob("~{in_output_prefix}/*")
        
    }
    runtime {
        cpu: in_threads
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/aardvark@sha256:aabe926c8022561fcf3f39f8c90e750335b1965897e1f517dc153c4eaf611a6b"
    }
}
