version 1.0

import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/vg_map_hts.wdl" as map
import "../tasks/vg_indexing.wdl" as index
import "../tasks/extract_reads.wdl" as extractReads_t

workflow HaplotypeSampling {
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
        description: "Create haplotype sampled graphs. It is a customized version of Parsa's wdl. It can take multiple haplotype numbers and create one graph per number. It does not create dist file."
    }

    parameter_meta {
        IN_GBZ_FILE: "Path to .gbz index file"
        INPUT_READ_FILE_ARRAY: "Input read files can be of any format like BAM, CRAM, FASTQ or FASTG.GZ"
        REFERENCE_FASTA: "Fasta file containing the reference genome for extracting reads (if for example the format is BAM)"
        HAPL_FILE: "Path to .hapl file"
        IN_DIST_FILE: "Path to .dist file"
        R_INDEX_FILE: "Path to .ri file"
        KFF_FILE: "Path to .kff file"
        IN_OUTPUT_NAME_PREFIX: "Name of the output file (Default: haplotype_sampled_graph)"
        IN_KMER_LENGTH: "Size of kmer using for sampling (Up to 31) (Default: 29)"
        CORES: "Number of cores to use with commands. (Default: 16)"
        WINDOW_LENGTH: "Window length used for building the minimizer index. (Default: 11)"
        SUBCHAIN_LENGTH: "Target length (in bp) for subchains. (Default: 10000)"
        HAPLOTYPE_NUMBER_ARRAY: "An array of numbers for generated synthetic haplotypes. (Default: [4])"
        PRESENT_DISCOUNT: "Multiplicative factor for discounting scores for present kmers. (Default: 0.9)"
        HET_ADJUST: "Additive term for adjusting scores for heterozygous kmers. (Default: 0.05)"
        ABSENT_SCORE: "Score for absent kmers. (Default: 0.8)"
        INCLUDE_REFERENCE: "Include reference paths and generic paths from the full graph in the sampled graph. (Default: true)"
        DIPLOID: "Activate diploid sampling. (Default: true)"
        CREATE_DISTANCE_INDEX_OPTIONS: ""
        DOCKER_IMAGE: "VG docker image"
    }
    input {
        File IN_GBZ_FILE
        Array[File] INPUT_READ_FILE_ARRAY
        File? REFERENCE_FASTA
        File? HAPL_FILE
        File? IN_DIST_FILE
        File? R_INDEX_FILE
        File? KFF_FILE
        String? IN_OUTPUT_NAME_PREFIX
        Int? IN_KMER_LENGTH
        Int CORES = 16
        Int WINDOW_LENGTH = 11
        Int SUBCHAIN_LENGTH = 10000
        Array[Int] HAPLOTYPE_NUMBER_ARRAY = [4]
        Float PRESENT_DISCOUNT = 0.9
        Float HET_ADJUST = 0.05
        Float ABSENT_SCORE = 0.8
        Boolean INCLUDE_REFERENCE = true
        Boolean DIPLOID = true
        String? SAMPLE_NAME_TO_REMOVE
        String DOCKER_IMAGE = "quay.io/vgteam/vg:v1.64.1"
    }

    String OUTPUT_NAME_PREFIX = select_first([IN_OUTPUT_NAME_PREFIX, "haplotype_sampled_graph"])
    Int KMER_LENGTH = select_first([IN_KMER_LENGTH, 29])

    # Have to create haplotype information
    if (!defined(HAPL_FILE)) {
        # create the dist index file and r-index file to create the haplotype information file .hapl

        if (!defined(IN_DIST_FILE)) {
            call index.createDistanceIndex {
                input:
                    in_gbz_file=IN_GBZ_FILE,
                    docker_image=DOCKER_IMAGE
            }
        }

        File dist_index_file = select_first([IN_DIST_FILE, createDistanceIndex.output_dist_index])

        if (!defined(R_INDEX_FILE)) {
            call index.createRIndex {
                input:
                    in_gbz_file=IN_GBZ_FILE,
                    nb_cores=CORES,
                    docker_image=DOCKER_IMAGE
            }
        }

        File r_index_file = select_first([R_INDEX_FILE, createRIndex.output_R_index])

        # create the haplotype information file
        call index.createHaplotypeIndex {
            input:
                in_gbz_file=IN_GBZ_FILE,
                in_dist_index=dist_index_file,
                in_R_index=r_index_file,
                nb_cores=CORES,
                kmer_length=KMER_LENGTH,
                window_length=WINDOW_LENGTH,
                subchain_length=SUBCHAIN_LENGTH,
                docker_image = DOCKER_IMAGE
        }
    }

    File haplotype_index = select_first([HAPL_FILE, createHaplotypeIndex.output_hap_index])

    # we make kmer database once and then 
    # run haplotype sampling for all haplotype numbers
    # using the same kmer database
    if (!defined(KFF_FILE)) {
        scatter (READ_FILE in INPUT_READ_FILE_ARRAY) {
            call extractReads_t.extractReads as extractReads {
                input:
                    readFile=READ_FILE,
                    referenceFasta=REFERENCE_FASTA,
                    memSizeGB=4,
                    threadCount=4,
                    diskSizeGB=floor(4 * size(READ_FILE, "GB")) + 32
            }
        }
        call utils.kmerCountingKMC {
            input:
                input_read_fastqs=extractReads.extractedRead,
                output_file_name=OUTPUT_NAME_PREFIX,
                kmer_length=KMER_LENGTH,
                nb_cores=CORES

        }

    }

    File kmer_information = select_first([KFF_FILE, kmerCountingKMC.kff_file])

    scatter (HAPLOTYPE_NUMBER in HAPLOTYPE_NUMBER_ARRAY) {
        call map.samplingHaplotypes {
            input:
                in_gbz_file=IN_GBZ_FILE,
                in_hap_index=haplotype_index,
                in_kmer_info=kmer_information,
                output_file_name = "${OUTPUT_NAME_PREFIX}.hap_num_${HAPLOTYPE_NUMBER}" ,
                haplotype_number = HAPLOTYPE_NUMBER,
                present_discount = PRESENT_DISCOUNT,
                het_adjust = HET_ADJUST,
                absent_score = ABSENT_SCORE,
                include_reference = INCLUDE_REFERENCE,
                nb_cores=CORES,
                use_diploid_sampling=DIPLOID,
                docker_image=DOCKER_IMAGE
        }
        if (defined(SAMPLE_NAME_TO_REMOVE)){
            call remove_sample_from_graph {
                input:
                    graph_gbz = samplingHaplotypes.output_graph,
                    sample_name = select_first([SAMPLE_NAME_TO_REMOVE]),
                    docker_image = DOCKER_IMAGE
            }
        }
        File sampled_graph_per_hap_number = select_first([remove_sample_from_graph.output_graph_gbz, samplingHaplotypes.output_graph])
    }


    output {
        Array[File] sampled_graph = sampled_graph_per_hap_number
    }

}


task remove_sample_from_graph {
    input {
        File graph_gbz
        String sample_name
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=128
        String docker_image="quay.io/vgteam/vg:v1.51.0"
        Int preemptible=2
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace
       
        INPUT_FILE=~{graph_gbz}
        INPUT_PREFIX=$(basename ${INPUT_FILE%%.gbz})
        OUTPUT_PREFIX="${INPUT_PREFIX}.~{sample_name}_removed"

        # create gbwt from input gbz
        vg gbwt \
            -Z ~{graph_gbz} \
            -o ${INPUT_PREFIX}.gbwt

        # remove sample from gbwt and create a new gbwt
        vg gbwt \
            ${INPUT_PREFIX}.gbwt \
            --remove-sample ~{sample_name}  \
            -o ${OUTPUT_PREFIX}.gbwt

        mkdir output
        # make a new gbz file using the new gbwt
        vg gbwt \
            ${OUTPUT_PREFIX}.gbwt \
            -x ~{graph_gbz} \
            --gbz-format  \
            --graph-name output/${OUTPUT_PREFIX}.gbz
 
    >>> 
    runtime {
        docker: docker_image
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File output_graph_gbz = glob("output/*.gbz")[0]
    }
}

