#!/bin/bash
qsub htseq-count -f bam -s no ./rna-seq_SRX263864.tophat2.bam ./homSap/Homo_sapiens.GRCh38.89.gtf.hg38 > ./rna-seq_SRX263864.tophat2.htseq.count
