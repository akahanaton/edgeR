#!/bin/bash
libSize=bowtie2.human.Aaron20170509.lib.size
if [ -e $libSize ]; then rm $libSize; fi
if [ -e libSize.human.sh ]; then rm libSize.human.sh; fi
#--------------------------------------------------
# for fq1 in `ls ./James_Transcriptome_Raw_Reads/*_R1.fastq`
#-------------------------------------------------- 
#--------------------------------------------------
# for fq1 in ./James_Transcriptome_Raw_Reads/JW_HEK293T*_R1.fastq.gz
#--------------------------------------------------
for fq1 in ./Aaron20170509/*[CW]*R1*.gz
do
    #--------------------------------------------------
    # ls $fq1
    #-------------------------------------------------- 
    fq2=`echo $fq1 | sed 's/_R1/_R2/'`
    id=`basename $fq1 _R1.fastq.gz`
    echo $id
    echo "#!/bin/bash" > $id.sh
    #--------------------------------------------------
    # echo "/gpfs/public/tools/bowtie2-2.2.5/bowtie2 --reorder -p 8 -x ./homSap/Homo_sapiens.GRCh38.dna_rm.toplevel.gtf.fa -1 $fq1 -2 $fq2 -S ./bowtie2/$id.sam > ./log/$id.bowtie2.log" >> $id.sh
    #--------------------------------------------------
    #--------------------------------------------------
    # echo "samtools view -b bowtie2/$id.sam -o ./bowtie2/$id.bam" >> $id.sh
    #--------------------------------------------------
    echo "/home/gmswenm/.local/bin/htseq-count -f bam -s no -o ./bowtie2/$id.bam.htseq -i gene_name ./bowtie2/$id.bam ./homSap/Homo_sapiens.GRCh38.84.gtf > ./bowtie2/$id.bam.count2" >> $id.sh
    if [[ $fq1 == *.gz ]]; then
        echo "zcat $fq1 | wc -l | awk -v id=$id '{print id,\$1/4}' >> $libSize" >> libSize.human.sh
    else
        echo "wc -l $fq1 | awk -v id=$id '{print id,\$1/4}' >> $libSize" >> libSize.human.sh
    fi
    chmod a+x $id.sh
    qsub -l h_vmem=40G -pe smp 8 ./$id.sh
done
chmod a+x libSize.human.sh
#--------------------------------------------------
# qsub bash ./libSize.human.sh
#--------------------------------------------------
