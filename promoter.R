library(GenomicRanges)
library(GenomicFeatures)
library(Biostrings)
library(Rsamtools)

#--------------------------------------------------
# txdb <- makeTxDbFromGFF("./pteAle/GCF_000325575.1_ASM32557v1_genomic.gff",format="gff3")
# seq.all <- readDNAStringSet("./pteAle/GCF_000325575.1_ASM32557v1_genomic.fa.mask/pteAle2.fa.masked", "fasta")
#--------------------------------------------------

upstream_pad<- 2000
downstream_pad<- 0

genes.all <- genes(txdb)

tss.all <- resize( genes.all, 1)

promoters.all <- promoters(tss.all, upstream_pad, downstream_pad)

genes.list = c("NFKB1","NFKB2","REL")
promoters.nfkb = promoters.all[promoters.all$gene_id %in% genes.list, ]

faFile = FaFile("./pteAle/GCF_000325575.1_ASM32557v1_genomic.fa.mask/pteAle2.fa.masked")

promoters.nfkb.seq<- getSeq(faFile, promoters.nfkb)

names(promoters.nfkb.seq) = genes.list

writeXStringSet(promoters.nfkb.seq, file="nfkb.fa")
