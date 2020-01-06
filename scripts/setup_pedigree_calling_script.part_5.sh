#!/bin/bash
#################################################################################################
##
##  Script to setup a bash script to run the 5th stage of the vg_wdl pedigree pipeline.
##
##  Inputs:
##
##  Assumptions:
##      The UDP cohort ID is the same as the UDP ID for the proband sample.
##
##  Last modified:
##  Last modified by: Charles Markello
##
#################################################################################################

## Create help statement
usage(){
cat << EOF

This script setups up a bash script to run a UDP cohort through the 5th stage of the vg_wdl
pedigree pipeline on the NIH Biowulf Cluster.

Inputs:
    -s List of Sibling UDP ID, Proband ID must be first in the list (in format UDP#### UDP#### UDP#### ...)
    -m Mother UDP ID (in format UDP####)
    -f Father UDP ID (in format UDP####)
    -w PATH to where the UDP cohort will be processed and where the input reads will be stored
    -g PATH to the workflow input directory
    -v PATH to the vg wdl repository
    -t (OPTIONAL, default=false) Set to 'true' if running workflow on small HG002 chr21 test data

Outputs:

Assumptions:

EOF

}

## Check number of arguments
if [ $# -lt 6 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## DEFAULT PARAMETERS
RUN_SMALL_TEST=false

## Parse through arguments
while getopts "s:m:f:w:g:v:t:h" OPTION; do
    case $OPTION in
        s)
            SIBLING_SAMPLE_NAMES+=($OPTARG)
        ;;
        m)
            MATERNAL_SAMPLE_NAME=$OPTARG
        ;;
        f)
            PATERNAL_SAMPLE_NAME=$OPTARG
        ;;
        w)
            COHORT_WORKFLOW_DIR=$OPTARG
        ;;
        g)
            WORKFLOW_INPUT_DIR=$OPTARG
        ;;
        v)
            VG_WDL_DIR=$OPTARG
        ;;
        t)
            RUN_SMALL_TEST=$OPTARG
        ;;
        h)
            usage
            exit 1
        ;;
        \?)
            usage
            exit 1
        ;;
    esac
done

PROBAND_SAMPLE_NAME="${SIBLING_SAMPLE_NAMES[0]}"
SIB_BAM_FILE_PARAMS=""
SIB_BAM_FILE_INDEX_PARAMS=""
SIB_ID_LIST=""
for SIBLING_ID in ${SIBLING_SAMPLE_NAMES[@]}
do
  SIB_BAM_FILE_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_sibling_map.final_outputs/output_links -name ${SIBLING_ID}_merged.positionsorted.bam))
  SIB_BAM_FILE_INDEX_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_sibling_map.final_outputs/output_links -name ${SIBLING_ID}_merged.positionsorted.bam.bai))
  SIB_BAM_FILE_PARAMS+="SIBLING_BAM_FILE_LIST='${SIB_BAM_FILE_PATH}' "
  SIB_BAM_FILE_INDEX_PARAMS+="SIBLING_BAM_FILE_INDEX_LIST='${SIB_BAM_FILE_INDEX_PATH}' "
  SIB_ID_LIST+="SAMPLE_NAME_SIBLING_LIST='${SIBLING_ID}' "
done
MATERNAL_GVCF_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_call.final_outputs/output_links -name ${MATERNAL_SAMPLE_NAME}_dragen_genotyped.gvcf.gz))
PATERNAL_GVCF_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_call.final_outputs/output_links -name ${PATERNAL_SAMPLE_NAME}_dragen_genotyped.gvcf.gz))
XG_FILE_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_parental_graph_construction.final_outputs/output_links -regex .*.xg))

rm -f ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_call.part_5.sh
echo '#!/bin/bash' >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_call.part_5.sh
echo "module load cromwell/40 python/3.6" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_call.part_5.sh
echo "source ${VG_WDL_DIR}/miniwdl_venv/bin/activate" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_call.part_5.sh
echo "cd ${COHORT_WORKFLOW_DIR}" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_call.part_5.sh
if [ $RUN_SMALL_TEST == false ]; then
    echo "miniwdl cromwell ${VG_WDL_DIR}/vg_wdl/workflows/vg_2nd_iter_pedigree_multi_call.wdl \\
MATERNAL_GVCF='${MATERNAL_GVCF_PATH}' \\
PATERNAL_GVCF='${PATERNAL_GVCF_PATH}' \\
${SIB_BAM_FILE_PARAMS} \\
${SIB_BAM_FILE_INDEX_PARAMS} \\
${SIB_ID_LIST} \\
PATH_LIST_FILE='${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt' \\
XG_FILE='${XG_FILE_PATH}' \\
REF_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa' \\
REF_INDEX_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai' \\
REF_DICT_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.dict' \\
VG_CONTAINER='quay.io/vgteam/vg:v1.19.0' \\
SNPEFF_DATABASE='${WORKFLOW_INPUT_DIR}/snpEff_v4_3_GRCh37.75.zip' \\
DRAGEN_REF_INDEX_NAME='hs37d5_v7' \\
UDPBINFO_PATH='Udpbinfo' \\
HELIX_USERNAME='${USER}' \\
DRAGEN_MODE='true' \\
SNPEFF_ANNOTATION='true' \\
-c ${VG_WDL_DIR}/vg_wdl/workflows/custom_biowulf_cromwell_singularity.conf \\
-d ${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_call.final_outputs" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_call.part_5.sh
else
    echo "miniwdl cromwell ${VG_WDL_DIR}/vg_wdl/workflows/vg_2nd_iter_pedigree_multi_call.wdl \\
MATERNAL_GVCF='${MATERNAL_GVCF_PATH}' \\
PATERNAL_GVCF='${PATERNAL_GVCF_PATH}' \\
${SIB_BAM_FILE_PARAMS} \\
${SIB_BAM_FILE_INDEX_PARAMS} \\
${SIB_ID_LIST} \\
PATH_LIST_FILE='${WORKFLOW_INPUT_DIR}/path_list_21.txt' \\
XG_FILE='${XG_FILE_PATH}' \\
REF_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa' \\
REF_INDEX_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai' \\
REF_DICT_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.dict' \\
VG_CONTAINER='quay.io/vgteam/vg:v1.19.0' \\
SNPEFF_DATABASE='${WORKFLOW_INPUT_DIR}/snpEff_v4_3_GRCh37.75.zip' \\
DRAGEN_REF_INDEX_NAME='hs37d5_v7' \\
UDPBINFO_PATH='Udpbinfo' \\
HELIX_USERNAME='${USER}' \\
DRAGEN_MODE='true' \\
SNPEFF_ANNOTATION='true' \\
-c ${VG_WDL_DIR}/vg_wdl/workflows/custom_biowulf_cromwell_singularity.conf \\
-d ${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_call.final_outputs" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_call.part_5.sh
fi

exit

