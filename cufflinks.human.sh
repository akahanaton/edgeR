#!/bin/bash
#--------------------------------------------------
# for fq1 in ./Aaron20170509/*[CW]*R1*.gz
#--------------------------------------------------
for fq1 in `ls -f ./James_Transcriptome_Raw_Reads/JW_HEK*_R1.fastq.gz`
do
    echo $fq1
    #--------------------------------------------------
    # ls $fq1
    #--------------------------------------------------
    fq2=`echo $fq1 | sed 's/R1/R2/'`
    id=`basename $fq1 _R1.fastq.gz`
    echo $id
    #--------------------------------------------------
    # ls $fq2
    #--------------------------------------------------
    qsub -pe smp 8 -N cuff.$id cufflinks -p 8 -b ./homSap/Homo_sapiens.GRCh38.dna_rm.toplevel.gtf.fa -G ./homSap/Homo_sapiens.GRCh38.84.gtf -g ./homSap/Homo_sapiens.GRCh38.84.gtf -o ./cufflinks/cufflinks.$id ./tophat2/tophat2.$id/accepted_hits.sam
done
