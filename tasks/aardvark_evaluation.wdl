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

task compareAardvarkSummaries {
    meta {
        description: "Join two Aardvark summary.tsv files into one table so the same metric from two runs can be read off the same line. The columns of summary.tsv are not hardcoded: any column whose values all parse as numbers is treated as a metric and gets its own left/right/delta columns, and the remaining columns are treated as the row's identity (variant type, stratification region, and so on) and used to match rows between the files."
    }
    input {
        File in_left_summary
        File in_right_summary
        String in_left_label = "left"
        String in_right_label = "right"
        String in_output_name = "comparison.tsv"
    }
    command <<<
        set -eux -o pipefail

        python3 <<'EOF'
        import csv

        def read_table(path):
            with open(path, newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                return list(reader.fieldnames or []), list(reader)

        left_columns, left_rows = read_table("~{in_left_summary}")
        right_columns, right_rows = read_table("~{in_right_summary}")

        # Preserve the left file's column order, and append any columns only the
        # right file has, so a change in Aardvark's output format still shows up.
        columns = left_columns + [c for c in right_columns if c not in left_columns]

        def is_number(text):
            try:
                float(text)
            except (TypeError, ValueError):
                return False
            return True

        # A column is a metric if nothing in it looks like a label. Empty cells
        # don't count against a column, so a metric that is missing for some
        # rows is still compared where it is present. The first column is always
        # kept as a label, so that rows can still be told apart in a table whose
        # row names happen to be numeric.
        def is_metric(column):
            if columns and column == columns[0]:
                return False
            values = [row.get(column) for row in left_rows + right_rows]
            values = [v for v in values if v not in (None, "")]
            return len(values) > 0 and all(is_number(v) for v in values)

        metric_columns = [c for c in columns if is_metric(c)]
        key_columns = [c for c in columns if c not in metric_columns]

        # With nothing to match rows on, fall back to position in the file.
        by_position = len(key_columns) == 0
        if by_position:
            key_columns = ["row"]

        def key_of(row, position):
            if by_position:
                return (str(position),)
            return tuple(row.get(c, "") or "" for c in key_columns)

        left_by_key = {key_of(row, i): row for i, row in enumerate(left_rows)}
        right_by_key = {key_of(row, i): row for i, row in enumerate(right_rows)}

        # Rows the two files don't share are kept rather than dropped; a run that
        # produced no calls at all in some stratum is exactly what we want to see.
        ordered_keys = list(left_by_key.keys())
        ordered_keys += [k for k in right_by_key.keys() if k not in left_by_key]

        header = list(key_columns)
        for column in metric_columns:
            header += ["~{in_left_label}_" + column, "~{in_right_label}_" + column, "delta_" + column]

        with open("~{in_output_name}", "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for key in ordered_keys:
                left_row = left_by_key.get(key, {})
                right_row = right_by_key.get(key, {})
                out = list(key)
                for column in metric_columns:
                    left_value = left_row.get(column, "") or ""
                    right_value = right_row.get(column, "") or ""
                    if is_number(left_value) and is_number(right_value):
                        delta = "%.10g" % (float(right_value) - float(left_value))
                    else:
                        delta = ""
                    out += [left_value, right_value, delta]
                writer.writerow(out)
        EOF
    >>>
    output {
        File output_comparison = "~{in_output_name}"
    }
    runtime {
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 10 SSD"
        docker: "python:3.12-slim"
    }
}