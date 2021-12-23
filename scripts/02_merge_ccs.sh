#!/bin/bash

# alorenzetti 202112
# description ####
# this script will merge partitioned
# ccs files generated separately

# requires pbbam from bioconda
# requires pbccs from bioconda

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate isoseqv3_env

for lib in ccfr ocfr ; do
    if [[ ! -d ./ccs_${lib} ]] ; then
        mkdir ./ccs_${lib}

        pbmerge -o ./ccs_${lib}/ccs_${lib}.bam ./ccs_12_parts_${lib}/${lib}_*.bam > ./ccs_${lib}/pbmerge.log 2> ./ccs_${lib}/pbmerge.err
        pbindex ./ccs_${lib}/ccs_${lib}.bam > ./ccs_${lib}/pbindex.log 2> ./ccs_${lib}/pbindex.err
    fi
done