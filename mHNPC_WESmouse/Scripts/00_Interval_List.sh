#!/bin/bash

################################################################################
#| Last update: 22/09/2025
#| Author: Ivana Rondon Lorefice
#|
#| Convert BED to Interval List Format
#|
#| Bed file was downloaded from: https://www.twistbioscience.com/resources/data-files/twist-mouse-exome-panel-bed-file
#| for Twist Mouse Exome Panel. This library is enriched following the annotations in the ensembl genome GRCm38. 
#| This bed file contains annotations of chromosomes starting with "chr"
#| I removed it before using the following command: 
#|       sed 's/^chr//' Twist_Mouse_Exome_Target_Rev1_7APR20.bed > Twist_Mouse_Exome_Target_Rev1_7APR20_Chr-Modified.bed
################################################################################

CONF_ENV="../00_conf_env.env"

source ${CONF_ENV}
source ${CONDA_ACTIVATE}
conda activate ${CONDA_ENV_WES}

picard BedToIntervalList I=${DATA_SHARED}${BED_FILE} \
    O=${DATA_SHARED}${INTERVAL_LIST} \
    SD=${GENOME_DICT}
