#!/bin/bash

# alorenzetti 202112

# description ####
# this script will take the
# ccs file and trim the adapters
# and set up the true orientation
# of transcripts using lima

# requires lima from bioconda

# setting up num threads
threads=15

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate isoseqv3_env

# requires lima from bioconda
for lib in ccfr ocfr ; do
    if [[ ! -d ./lima_${lib} ]] ; then
        mkdir ./lima_${lib}

        # observing the ccs reads I was able to detect
        # the adaps are those from NEB reported in
        # https://github.com/PacificBiosciences/IsoSeq/blob/master/isoseq-clustering.md
        # I will create the required files
        echo -e ">NEB_5p\nGCAATGAAGTCGCAGGGTTGGGG\n>NEB_3p\nGTACTCTGCGTTGATACCACTGCTT" > ./lima_${lib}/adaps.fa

        # performing lima
        lima ./ccs_${lib}/ccs_${lib}.bam \
            ./lima_${lib}/adaps.fa \
            ./lima_${lib}/ccs_${lib}_lima.bam \
            --isoseq --peek-guess -j $threads > ./lima_${lib}/lima.log 2> ./lima_${lib}/lima.err
    fi
done