#!/bin/bash

# alorenzetti 202112

# description ####
# this script will take the
# ccs trimmed files and remove
# poly-a and concatamers by
# using isoseq refine bringing up
# flnc tx (full-length non-concatamer)
# then output files
# will be submitted to 
# the clustering step

# requires isoseq3 from bioconda

# setting up num threads
threads=15

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate isoseqv3_env

for lib in ccfr ocfr ; do
    if [[ ! -d ./isoseq_refine_cluster_${lib} ]] ; then
        mkdir ./isoseq_refine_cluster_${lib}

        isoseq3 refine -j $threads --require-polya \
                    ./lima_${lib}/ccs_${lib}_lima.NEB_5p--NEB_3p.bam \
                    ./lima_${lib}/adaps.fa \
                    ./isoseq_refine_cluster_${lib}/ccs_${lib}_lima_refine_onlyPolyA.bam > ./isoseq_refine_cluster_${lib}/refine_onlyPolyA.log 2> ./isoseq_refine_cluster_${lib}/refine_onlyPolyA.err

        isoseq3 cluster -j $threads --use-qvs --verbose \
                        ./isoseq_refine_cluster_${lib}/ccs_${lib}_lima_refine_onlyPolyA.bam \
                        ./isoseq_refine_cluster_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster.bam > ./isoseq_refine_cluster_${lib}/cluster_onlyPolyA.log 2> ./isoseq_refine_cluster_${lib}/cluster_onlyPolyA.err

    fi
done