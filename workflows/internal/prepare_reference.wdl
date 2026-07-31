version 1.0

import "../../tasks/bioinfo_utils.wdl" as utils
import "../../tasks/validation.wdl" as validation
import "../../tasks/vg_map_hts.wdl" as map

workflow PrepareReference {
    meta {
        description: "## Prepare reference workflow \n Work out which contigs to work on and get a FASTA reference for them, with its .fai and .dict. Anything the caller provides is used as it is; a contig list and a reference that aren't provided are extracted from the graph, so that the reference's contigs match the graph's paths exactly. This is what the mapping and calling workflows need before they can surject, split by contig, or call variants; it is a workflow of its own so that a caller running several pipelines over one graph can do it once and keep every run on the same contigs and the same reference bases."
    }

    parameter_meta {
        GBZ_FILE: "Path to the .gbz index file to take contigs and reference bases from"
        CONTIGS: "(OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index. If using REFERENCE_PREFIX, contig names in here should have the prefix."
        PATH_LIST_FILE: "(OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths. If using REFERENCE_PREFIX, contig names in here should have the prefix."
        REFERENCE_PREFIX: "Prefix to remove from path names to get reference contig names. If not empty, an extracted contig list is checked to make sure it carries the prefix."
        REFERENCE_FILE: "(OPTIONAL) If specified, use this FASTA reference, gzipped or not, instead of extracting it from the graph. Required if the graph does not contain all bases of the reference. If using REFERENCE_PREFIX, contig names in here should not have the prefix."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths."
        EXTRACT_MEM: "Memory, in GB, to use when extracting contigs and reference bases from the graph. Default is 120."
        VG_DOCKER: "Container image to use when running vg"
    }

    input {
        File GBZ_FILE
        Array[String]+? CONTIGS
        File? PATH_LIST_FILE
        String REFERENCE_PREFIX = ""
        File? REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        File? REFERENCE_DICT_FILE
        Int EXTRACT_MEM = 120
        String VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
    }

    if (!defined(CONTIGS)) {
        if (!defined(PATH_LIST_FILE)) {
            # Extract path names to call against from the GBZ file if neither a
            # path list nor a contig list was provided.
            # Filter down to major paths, because GRCh38 includes thousands of
            # decoys and unplaced/unlocalized contigs, and we can't efficiently
            # scatter across them, nor do we care about accuracy on them, and also
            # calling on the decoys is semantically meaningless.
            call map.extractSubsetPathNames {
                input:
                    in_gbz_file=GBZ_FILE,
                    in_reference_prefix=REFERENCE_PREFIX,
                    in_extract_mem=EXTRACT_MEM,
                    vg_docker=VG_DOCKER
            }

            if (REFERENCE_PREFIX != "") {
                # A list we made ourselves can still come out wrong, if the
                # graph's paths aren't named the way the prefix says they are.
                call validation.checkPathList as checkExtractedPathList {
                    input:
                        in_path_list_file=extractSubsetPathNames.output_path_list_file,
                        in_reference_prefix=REFERENCE_PREFIX
                }
            }
        }
    }
    if (defined(CONTIGS)) {
        # Put the paths in a file to use later. We know the value is defined,
        # but WDL is a bit low on unboxing calls for optionals so we use
        # select_first.
        File written_path_names_file = write_lines(select_first([CONTIGS]))
    }
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, extractSubsetPathNames.output_path_list_file, written_path_names_file])

    # To make sure that we have a FASTA reference with a contig set that
    # exactly matches the graph (except for removing the name prefix), we
    # generate it ourselves, from the graph.
    if (!defined(REFERENCE_FILE)) {
        call map.extractReference {
            input:
            in_gbz_file=GBZ_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_prefix_to_strip=REFERENCE_PREFIX,
            in_extract_mem=EXTRACT_MEM,
            vg_docker=VG_DOCKER
        }
    }
    if (defined(REFERENCE_FILE)) {
        call utils.uncompressReferenceIfNeeded {
            input:
            # We know REFERENCE_FILE is defined but the WDL type system doesn't.
            in_reference_file=select_first([REFERENCE_FILE])
        }
    }
    File resolved_reference_file = select_first([uncompressReferenceIfNeeded.reference_file, extractReference.reference_file])

    # The .fai and the .dict are made by the same task, so a caller that brought
    # only one of them still needs it run.
    if (!defined(REFERENCE_INDEX_FILE) || !defined(REFERENCE_DICT_FILE)) {
        call utils.indexReference {
            input:
                in_reference_file=resolved_reference_file
        }
    }

    output {
        File path_list_file = pipeline_path_list_file
        File reference_file = resolved_reference_file
        File reference_index_file = select_first([REFERENCE_INDEX_FILE, indexReference.reference_index_file])
        File reference_dict_file = select_first([REFERENCE_DICT_FILE, indexReference.reference_dict_file])
    }
}
