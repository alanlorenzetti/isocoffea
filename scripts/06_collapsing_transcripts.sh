#!/bin/bash

# alorenzetti 202112

# description ####
# this script will run cDNA cupcake
# in order to collapse clusters
# output by isoseq3

# requires cDNA cupcake and its dependencies
# cDNA cupcake should be installed following 
# https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake:-supporting-scripts-for-Iso-Seq-after-clustering-step#install

# requires samtools

# setting up num threads
threads=15

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate cogent_env

for lib in ccfr ocfr ; do
    if [[ ! -d collapse_isoforms_${lib} ]] ; then
        mkdir collapse_isoforms_${lib}

        # downloading ref genome and annotation
        # from NCBI assembly
        # https://www.ncbi.nlm.nih.gov/assembly/GCF_003713225.1
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_genomic.fna.gz -O collapse_isoforms_${lib}/Carabica.fa.gz > /dev/null 2>&1
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_genomic.gff.gz -O collapse_isoforms_${lib}/Carabica.gff.gz > /dev/null 2>&1

        # gunzipping
        gunzip collapse_isoforms_${lib}/Carabica.fa.gz
        gunzip collapse_isoforms_${lib}/Carabica.gff.gz

        # copying clustered transcripts file to target directory
        cp isoseq_refine_cluster_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster.hq.fasta.gz collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster.hq.fasta.gz
        gunzip collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster.hq.fasta.gz

        # aligning reads to the reference genome
        minimap2 -a \
                 -x splice \
                 -t $threads \
                 -u f \
                 --secondary=no \
                 -C5 \
                 collapse_isoforms_${lib}/Carabica.fa \
                 collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster.hq.fasta > collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq.sam

        sort -k 3,3 -k 4,4n collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq.sam > collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_sorted.sam

        # creating a version with supp alignments for visualization purposes
        minimap2 -a \
                 -x splice \
                 -t $threads \
                 -u f \
                 -C5 \
                 collapse_isoforms_${lib}/Carabica.fa \
                 collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster.hq.fasta > collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_viz.sam

        samtools sort -@ $threads \
                 -o collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_viz.bam \
                 collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_viz.sam

        samtools index -b collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_viz.bam

        # collapsing
        collapse_isoforms_by_sam.py --input collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster.hq.fasta \
                                    -s collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_sorted.sam \
                                    --dun-merge-5-shorter \
                                    -o collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_sorted \
                                    > collapse_isoforms_${lib}/collapse_isoforms_by_sam_out.log \
                                    2> collapse_isoforms_${lib}/collapse_isoforms_by_sam_err.log

        # getting associated abundance
        get_abundance_post_collapse.py collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_sorted.collapsed isoseq_refine_cluster_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster.cluster_report.csv

        # filter away 5' degraded isoforms
        filter_away_subset.py collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_sorted.collapsed

        # aligning again to check results
        # using a genome browser
        minimap2 -a \
                 -x splice \
                 -t $threads \
                 -u f \
                 -C5 \
                 collapse_isoforms_${lib}/Carabica.fa \
                 collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_sorted.collapsed.rep.fa > collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_sorted_collapsed_filtered_rep.sam

        samtools sort -@ $threads \
                 -o collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_sorted_collapsed_filtered_rep.bam \
                 collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_sorted_collapsed_filtered_rep.sam

        samtools index -b collapse_isoforms_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_sorted_collapsed_filtered_rep.bam
    fi
done