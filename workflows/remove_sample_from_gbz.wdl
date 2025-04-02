version 1.0

import "haplotype_sampling_customized.wdl" as hap_wdl
import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/vg_indexing.wdl" as index

workflow RemoveSampleFromGraph {
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
        description: "Remove a sample from gbz index file and create dist, r and haplotype indices for the output gbz"
    }

    parameter_meta {
        IN_GBZ_FILE: "Path to .gbz index file"
        SAMPLE_NAME_TO_REMOVE: "Name of the sample to remove from graph"
        IN_OUTPUT_NAME_PREFIX: "Name of the output file (Default: haplotype_sampled_graph)"
        IN_KMER_LENGTH: "Size of kmer using for sampling (Up to 31) (Default: 29)"
        CORES: "Number of cores to use with commands. (Default: 16)"
        WINDOW_LENGTH: "Window length used for building the minimizer index. (Default: 11)"
        SUBCHAIN_LENGTH: "Target length (in bp) for subchains. (Default: 10000)"
        DOCKER_IMAGE: "VG docker image"
    }
    input {
        File IN_GBZ_FILE
        String SAMPLE_NAME_TO_REMOVE
        String? IN_OUTPUT_NAME_PREFIX
        Int? IN_KMER_LENGTH
        Int CORES = 16
        Int WINDOW_LENGTH = 11
        Int SUBCHAIN_LENGTH = 10000
        String DOCKER_IMAGE = "quay.io/vgteam/vg:v1.64.1"
    }

    String OUTPUT_NAME_PREFIX = select_first([IN_OUTPUT_NAME_PREFIX, "haplotype_sampled_graph"])
    Int KMER_LENGTH = select_first([IN_KMER_LENGTH, 29])

    call hap_wdl.remove_sample_from_graph {
        input:
              graph_gbz = IN_GBZ_FILE,
              sample_name = select_first([SAMPLE_NAME_TO_REMOVE]),
              docker_image = DOCKER_IMAGE
    }

    call index.createDistanceIndex {
        input:
            in_gbz_file = remove_sample_from_graph.output_graph_gbz,
            docker_image = DOCKER_IMAGE
    }
    
    call index.createRIndex {
        input:
            in_gbz_file = remove_sample_from_graph.output_graph_gbz,
            nb_cores = CORES,
            docker_image = DOCKER_IMAGE
    }

    # create the haplotype information file
    call index.createHaplotypeIndex {
        input:
            in_gbz_file = remove_sample_from_graph.output_graph_gbz,
            in_dist_index = createDistanceIndex.output_dist_index,
            in_R_index = createRIndex.output_R_index,
            nb_cores = CORES,
            kmer_length = KMER_LENGTH,
            window_length = WINDOW_LENGTH,
            subchain_length = SUBCHAIN_LENGTH,
            docker_image = DOCKER_IMAGE
     }

     output {
         File output_gbz = remove_sample_from_graph.output_graph_gbz
         File dist_index = createDistanceIndex.output_dist_index
         File r_index = createRIndex.output_R_index
         File haplotype_index = createHaplotypeIndex.output_hap_index 
     }
}
