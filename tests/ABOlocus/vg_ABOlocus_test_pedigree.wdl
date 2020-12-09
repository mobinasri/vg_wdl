version 1.0

# This example uses vg_construct_and_index.wdl to build a graph including the
# ABO locus from GRCh38 (a 50Kbp region including the gene), with the 1000
# Genomes small variants (SNVs & short indels) and individual haplotypes.
# It then runs through and tests the vg pedigree workflow
# using vg_trio_multi_map_call.wdl against HG002 trio reads.

import "../../workflows/vg_construct_and_index.wdl"
import "../../workflows/vg_trio_multi_map_call.wdl"

workflow vg_ABOlocus_test {
    input {
        File ABOlocus_fa_gz
        File ABOlocus_small_vcf_gz
        File maternal_reads_bam
        File paternal_reads_bam
        File child_reads_bam
        File ref_file
        File ref_index_file
        File ref_dict_file
        File ref_file_gz
        File ped_file
        String vg_docker = 'quay.io/vgteam/vg:v1.28.0'
    }

    # build & check the ABOlocus graph
    call vg_construct_and_index.vg_construct_and_index as cons { input:
        graph_name = "ABOlocus_small",
        ref_fasta_gz = ABOlocus_fa_gz,
        contigs = ["ABOlocus"],
        contigs_vcf_gz = [ABOlocus_small_vcf_gz],
        giraffe_indexes = true,
        vg_docker = vg_docker
    }

    # extract FASTQs from trio BAMs for use with vg trio workflow
    call bam_to_paired_fastq as child_reads { input:
        bam = child_reads_bam
    }
    call bam_to_paired_fastq as mom_reads { input:
        bam = maternal_reads_bam
    }
    call bam_to_paired_fastq as dad_reads { input:
        bam = paternal_reads_bam
    }
    
    call vg_trio_multi_map_call.vgTrioPipeline as vgtrio { input:
        MATERNAL_INPUT_READ_FILE_1 = mom_reads.fastq_1_gz,
        MATERNAL_INPUT_READ_FILE_2 = mom_reads.fastq_2_gz,
        PATERNAL_INPUT_READ_FILE_1 = dad_reads.fastq_1_gz,
        PATERNAL_INPUT_READ_FILE_2 = dad_reads.fastq_2_gz,
        SIBLING_INPUT_READ_FILE_1_LIST = [child_reads.fastq_1_gz],
        SIBLING_INPUT_READ_FILE_2_LIST = [child_reads.fastq_2_gz],
        SAMPLE_NAME_SIBLING_LIST = ["HG002"],
        SAMPLE_NAME_MATERNAL = "HG004",
        SAMPLE_NAME_PATERNAL = "HG003",
        XG_FILE = cons.xg,
        GBWT_FILE = cons.gbwt,
        GGBWT_FILE = cons.ggbwt,
        DIST_FILE = cons.dist,
        MIN_FILE = cons.min,
        GRAPH_NAME = "ABOlocus_parental",
        CONTIGS = ["ABOlocus"],
        REF_FILE = ref_file,
        REF_INDEX_FILE = ref_index_file,
        REF_DICT_FILE = ref_dict_file,
        SPLIT_READ_CORES = 1,
        SPLIT_READ_DISK = 5,
        MAP_CORES = 1,
        MAP_DISK = 5,
        MAP_MEM = 4,
        VGCALL_CORES = 1,
        VGCALL_DISK = 5,
        VGCALL_MEM = 4,
        REF_FASTA_GZ = ref_file_gz,
        PED_FILE = ped_file,
        SNPEFF_ANNOTATION = false,
        CLEANUP_FILES = false,
        USE_DECOYS = false,
        MAPPER="GIRAFFE",
        VG_CONTAINER = vg_docker 
    }
    
    call check_trio_bams { input:
        mom_fastq_1_gz = mom_reads.fastq_1_gz,
        dad_fastq_1_gz = dad_reads.fastq_1_gz,
        child_fastq_1_gz = child_reads.fastq_1_gz,
        mom_bam = vgtrio.output_maternal_bam,
        dad_bam = vgtrio.output_paternal_bam,
        child_bam = vgtrio.final_output_sibling_bam_list[0]
    }
    
    output {
        File output_cohort_vcf = vgtrio.output_cohort_vcf
        Int n_mom_reads = check_trio_bams.n_mom_reads
        Int n_dad_reads = check_trio_bams.n_dad_reads
        Int n_child_reads = check_trio_bams.n_child_reads
        Int n_mom_aligned_reads = check_trio_bams.n_mom_aligned_reads
        Int n_dad_aligned_reads = check_trio_bams.n_dad_aligned_reads
        Int n_child_aligned_reads = check_trio_bams.n_child_aligned_reads
        Int mom_reads_aligned_identically = check_trio_bams.mom_reads_aligned_identically
        Int dad_reads_aligned_identically = check_trio_bams.dad_reads_aligned_identically
        Int child_reads_aligned_identically = check_trio_bams.child_reads_aligned_identically
    }
}

task bam_to_paired_fastq {
    input {
        File bam
    }

    command <<<
        set -eux -o pipefail
        
        samtools sort -n -@ $(nproc) -O BAM -o namesorted.bam "~{bam}"
        nm=$(basename "~{bam}" .bam)
        java -jar /picard.jar SamToFastq I=namesorted.bam RE_REVERSE=true INCLUDE_NON_PF_READS=true "FASTQ=${nm}_1.fastq" "SECOND_END_FASTQ=${nm}_2.fastq" "UNPAIRED_FASTQ=${nm}_unpaired.fastq" VALIDATION_STRINGENCY=LENIENT
        pigz *.fastq
    >>>

    runtime {
        docker: "quay.io/cmarkello/bamsplit@sha256:5ebf24ab2647f481cdb2e0827aea12ba60ae7ede9eda9cb97ce58b5e6b7e2e0a"
    }

    output {
        File fastq_1_gz = glob("*_1.fastq.gz")[0]
        File fastq_2_gz = glob("*_2.fastq.gz")[0]
        File fastq_unpaired_gz = glob("*_unpaired.fastq.gz")[0]
    }
}

task check_trio_bams {
    # checks if the trio bams have the expected # of mappings, and finds the proportion aligned at 100% identity
    input {
        File mom_fastq_1_gz
        File dad_fastq_1_gz
        File child_fastq_1_gz
        File? mom_bam
        File? dad_bam
        File? child_bam
    }

    command <<<
        set -eux -o pipefail
        
        n_mom_reads=$(expr $(zcat ~{mom_fastq_1_gz} | wc -l) / 2)
        n_mom_aligned_reads=$(samtools view ~{mom_bam} | wc -l)
        diff_mom_aligned_source_reads=$(($n_mom_reads - $n_mom_aligned_reads))
        if [ ${diff_mom_aligned_source_reads#-} -gt 10 ]; then
            echo "wrong read count for maternal alignments" >&2
            exit 1
        fi
        echo "$n_mom_reads" > n_mom_reads
        echo "$n_mom_aligned_reads" > n_mom_aligned_reads
        samtools view "~{mom_bam}" | perl -lane 'print if $F[5] =~ /^250M$/;' | wc -l > n_mom_identical
        
        n_dad_reads=$(expr $(zcat ~{dad_fastq_1_gz} | wc -l) / 2)
        n_dad_aligned_reads=$(samtools view ~{dad_bam} | wc -l)
        diff_dad_aligned_source_reads=$(($n_dad_reads - $n_dad_aligned_reads))
        if [ ${diff_dad_aligned_source_reads#-} -gt 10 ]; then
            echo "wrong read count for paternal alignments" >&2
            exit 1
        fi
        echo "$n_dad_reads" > n_dad_reads
        echo "$n_dad_aligned_reads" > n_dad_aligned_reads
        samtools view "~{dad_bam}" | perl -lane 'print if $F[5] =~ /^250M$/;' | wc -l > n_dad_identical
        
        n_child_reads=$(expr $(zcat ~{child_fastq_1_gz} | wc -l) / 2)
        n_child_aligned_reads=$(samtools view ~{child_bam} | wc -l)
        diff_child_aligned_source_reads=$(($n_child_reads - $n_child_aligned_reads))
        if [ ${diff_child_aligned_source_reads#-} -gt 15 ]; then
            echo "wrong read count for child alignments" >&2
            exit 1
        fi
        echo "$n_child_reads" > n_child_reads
        echo "$n_child_aligned_reads" > n_child_aligned_reads
        samtools view "~{child_bam}" | perl -lane 'print if $F[5] =~ /^250M$/;' | wc -l > n_child_identical
    >>>

    runtime {
        docker: "quay.io/cmarkello/bamsplit@sha256:5ebf24ab2647f481cdb2e0827aea12ba60ae7ede9eda9cb97ce58b5e6b7e2e0a"
    }

    output {
        Int n_mom_reads = read_int("n_mom_reads")
        Int n_dad_reads = read_int("n_dad_reads")
        Int n_child_reads = read_int("n_child_reads")
        Int n_mom_aligned_reads = read_int("n_mom_aligned_reads")
        Int n_dad_aligned_reads = read_int("n_dad_aligned_reads")
        Int n_child_aligned_reads = read_int("n_child_aligned_reads")
        Int mom_reads_aligned_identically = read_int("n_mom_identical")
        Int dad_reads_aligned_identically = read_int("n_dad_identical")
        Int child_reads_aligned_identically = read_int("n_child_identical")
    }
}

