#!/bin/bash

# alorenzetti 202112
# description ####
# this script will perform the isoseq v3 pipeline
# according to
# https://github.com/PacificBiosciences/IsoSeq/blob/master/isoseq-clustering.md

# requires pbccs from bioconda

# usage 
# bash process_isoseq.sh yes 12

# we need to run the 12/12 parts

# setting vars ####
threads=15

# getting args
# use chunking
chunking=$1

# total chunks
chunk=$2

# total chunks
totchunks=$chunk

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate isoseqv3_env

# file names for CCFR and OCFR libs
ccfr=/home/alorenzetti/quillaja_bucket/coffea_isoseq/raw/m64165_211109_030649.subreads.bc1013--bc1013.bam
ocfr=/home/alorenzetti/quillaja_bucket/coffea_isoseq_OCFR/raw/m64165_211109_030649.subreads.bc1014--bc1014.bam

# creating directory to store results
if [[ $chunking == "yes" ]] ; then
    if [[ ! -d ./ccs_${totchunks}_parts_ccfr ]] ; then mkdir ./ccs_${totchunks}_parts_ccfr ; fi
    if [[ ! -d ./ccs_${totchunks}_parts_ocfr ]] ; then mkdir ./ccs_${totchunks}_parts_ocfr ; fi

    for i in $(seq 1 ${totchunks}) ; do
        if [[ ! -f ./ccs_${totchunks}_parts_ccfr/ccfr_${i}.bam ]] ; then
            ccs \
            $ccfr \
            ./ccs_${totchunks}_parts_ccfr/ccfr_${i}.bam \
            --chunk ${i}/${totchunks} -j $threads 
        fi

        if [[ ! -f ./ccs_${totchunks}_parts_ocfr/ocfr_${i}.bam ]] ; then
            ccs \
            $ocfr \
            ./ccs_${totchunks}_parts_ocfr/ocfr_${i}.bam \
            --chunk ${i}/${totchunks} -j $threads 
        fi
    done

else
    echo "The program does not support $chunking as a valid argument anymore."
    exit 1
fi