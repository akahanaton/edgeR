#!/bin/bash
for fq1 in ./James_Transcriptome_Raw_Reads/JW_HEK293T*R1.fastq.gz
do
    #--------------------------------------------------
    # ls $fq1
    #-------------------------------------------------- 
    fq2=`echo $fq1 | sed 's/R1/R2/'`
    #--------------------------------------------------
    # ls $fq2
    #-------------------------------------------------- 
    id=`basename $fq1 _R1.fastq.gz`
    qsub -pe smp 8 tophat2 -p 8 --no-convert-bam -G ./homSap/Homo_sapiens.GRCh38.84.gtf -r 300 --mate-std-dev 100 -o tophat2/tophat2.$id ./homSap/Homo_sapiens.GRCh38.dna_rm.toplevel.gtf.fa $fq1 $fq2
    #--------------------------------------------------
    # /gpfs/public/tools2/cufflinks-2.2.1.Linux_x86_64/cufflinks -p 25 -G ./Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.84.gtf -o cufflinks.$id ./tophat2.$id/accepted_hits.bam
    #--------------------------------------------------
done
