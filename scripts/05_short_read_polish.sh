#!/bin/bash

# alorenzetti 202112

# description ####
# this script will take
# the FLNC clustered reads
# and polish using a short read
# dataset obtained from similar samples

# requires trimmomatic from bioconda: trimmomatic_env
# requires lordec from bioconda: lordec_env
# requires seqtk from bioconda: seqtk_env
# requires pigz: seqtk_env

# setting up num threads
threads=15

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
# activating trimmomatic env
conda activate trimmomatic_env

for lib in ccfr ocfr ; do
    if [[ $lib == "ccfr" ]] ; then libalt=CC ; else libalt=OC ; fi

    if [[ ! -d short_read_polish_${lib} ]] ; then
        mkdir short_read_polish_${lib}

        # creating adapter fasta according to
        # facility report
        echo -e ">read1_adap\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCA" > short_read_polish_${lib}/adap.fa
        echo -e ">read2_adap\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" >> short_read_polish_${lib}/adap.fa
    fi

    if [[ ! -d short_read_polish_${lib}/trimmed ]] ; then
        mkdir short_read_polish_${lib}/trimmed

        # the first step is to trim short-read
        # dataset to include only high quality reads

        # running trimmomatic for samples
        prefixes=`ls /home/alorenzetti/quillaja_bucket/coffea_hx/raw/LCS7609_DS_${libalt}_leaf*.fq.gz | sed 's/_R[12].fq.gz//' | sort | uniq`

        for prefix in $prefixes ; do
            sample=${prefix##*/}

            trimmomatic PE \
                        -threads $threads \
                        ${prefix}_R1.fq.gz ${prefix}_R2.fq.gz \
                        short_read_polish_${lib}/trimmed/${sample}-paired_1.fq.gz short_read_polish_${lib}/trimmed/${sample}-unpaired_1.fq.gz \
                        short_read_polish_${lib}/trimmed/${sample}-paired_2.fq.gz short_read_polish_${lib}/trimmed/${sample}-unpaired_2.fq.gz \
                        ILLUMINACLIP:short_read_polish_${lib}/adap.fa:1:30:10 \
                        SLIDINGWINDOW:4:30 \
                        MINLEN:50 > short_read_polish_${lib}/trimmed/${sample}_out.log 2> short_read_polish_${lib}/trimmed/${sample}_err.log
        done
    fi

    # activating seqtk env
    # conda activate seqtk_env

    # # subsetting short reads is not necessary for this case
    # # we cannot use more than 128GB RAM
    # # using 1/10 of the reads
    # prefixes=`ls short_read_polish_${lib}/trimmed/*_${libalt}*.fq.gz | sed 's/-.*_[12].fq.gz//' | sort | uniq`
    # for prefix in $prefixes ; do 
    #     if [[ -d short_read_polish_${lib}/trimmed ]] ; then 
    #         if [[ ! -f ${prefix}-paired-sub_1.fq.gz || ! -f ${prefix}-paired-sub_2.fq.gz ]] ; then 
    #             lines=`unpigz -p $threads -c ${prefix}-paired_1.fq.gz | wc -l`
    #             subsample=$((lines/4/4))

    #             seqtk sample -s 665 ${prefix}-paired_1.fq.gz $subsample | pigz -p $threads -c > ${prefix}-paired-sub_1.fq.gz
    #             seqtk sample -s 665 ${prefix}-paired_2.fq.gz $subsample | pigz -p $threads -c > ${prefix}-paired-sub_2.fq.gz
    #         fi

    #         if [[ ! -f ${prefix}-unpaired-sub_1.fq.gz ]] ; then 
    #             lines=`unpigz -p $threads -c ${prefix}-unpaired_1.fq.gz | wc -l`
    #             subsample=$((lines/4/4))

    #             seqtk sample -s 665 ${prefix}-unpaired_1.fq.gz $subsample | pigz -p $threads -c > ${prefix}-unpaired-sub_1.fq.gz
    #         fi

    #         if [[ ! -f ${prefix}-unpaired-sub_2.fq.gz ]] ; then 
    #             lines=`unpigz -p $threads -c ${prefix}-unpaired_2.fq.gz | wc -l`
    #             subsample=$((lines/4/4))

    #             seqtk sample -s 665 ${prefix}-unpaired_2.fq.gz $subsample | pigz -p $threads -c > ${prefix}-unpaired-sub_2.fq.gz
    #         fi
    #     fi
    # done

    # activating lordec env
    conda activate lordec_env

    # the next step is to run fmlrc to correct
    # clustered FLNC reads
    if [[ ! -d short_read_polish_${lib}/polished ]] ; then
        mkdir short_read_polish_${lib}/polished

        cd short_read_polish_${lib}/polished

        trimmedreads=`ls ../trimmed/*${libalt}*[12].fq.gz | xargs | sed 's/ /,/g'`

        lordec-correct -T $threads \
                    -2 $trimmedreads \
                    -k 15 \
                    -s 3 \
                    -i /home/alorenzetti/isoseq_refine_cluster_${lib}/ccs_${lib}_lima_refine_onlyPolyA_cluster.hq.fasta.gz \
                    -o ccs_${lib}_lima_refine_onlyPolyA_cluster_hq_polished.fasta \
                    > lordec_out.log 2> lordec_err.log

        cd ../..
    fi

done
