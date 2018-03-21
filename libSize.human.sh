zcat ./Aaron20170509/Y10-W138PD26_S10_L001_R1_001.fastq.gz | wc -l | awk -v id=Y10-W138PD26_S10_L001 '{print id,$1/4}' >> bowtie2.human.Aaron20170509.lib.size
zcat ./Aaron20170509/Y11-W138PD30_S11_L001_R1_001.fastq.gz | wc -l | awk -v id=Y11-W138PD30_S11_L001 '{print id,$1/4}' >> bowtie2.human.Aaron20170509.lib.size
zcat ./Aaron20170509/Y12-W138PD35_S12_L001_R1_001.fastq.gz | wc -l | awk -v id=Y12-W138PD35_S12_L001 '{print id,$1/4}' >> bowtie2.human.Aaron20170509.lib.size
zcat ./Aaron20170509/Y7-CR560495_S7_L001_R1_001.fastq.gz | wc -l | awk -v id=Y7-CR560495_S7_L001 '{print id,$1/4}' >> bowtie2.human.Aaron20170509.lib.size
zcat ./Aaron20170509/Y8-CR560947_S8_L001_R1_001.fastq.gz | wc -l | awk -v id=Y8-CR560947_S8_L001 '{print id,$1/4}' >> bowtie2.human.Aaron20170509.lib.size
zcat ./Aaron20170509/Y9-CR560100_S9_L001_R1_001.fastq.gz | wc -l | awk -v id=Y9-CR560100_S9_L001 '{print id,$1/4}' >> bowtie2.human.Aaron20170509.lib.size
