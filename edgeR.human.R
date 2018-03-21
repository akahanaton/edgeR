library(edgeR)
human0 <- read.table("./bowtie2/JW_HEK293T_Mock_0h.sam.count2", header=F,sep="\t", stringsAsFactors=F )
human8 <- read.table("./bowtie2/JW_HEK293T_HeV_8h.sam.count2", header=F,sep="\t", stringsAsFactors=F )
human24 <- read.table("./bowtie2/JW_HEK293T_HeV_24h.sam.count2", header=F,sep="\t", stringsAsFactors=F )
human0_total <- read.table("./James_Transcriptome_Raw_Reads/JW_HEK293T_Mock_0h_R1.fastq.readcount", header=F,sep=" ", stringsAsFactors=F)
human8_total <- read.table("./James_Transcriptome_Raw_Reads/JW_HEK293T_HeV_8h_R1.fastq.readcount", header=F,sep=" ", stringsAsFactors=F)
human24_total <- read.table("./James_Transcriptome_Raw_Reads/JW_HEK293T_HeV_24h_R1.fastq.readcount", header=F,sep=" ", stringsAsFactors=F)

human_total = c(human0_total[,3], human8_total[,3], human24_total[,3])
human = cbind(human0, human8[,2], human24[,2])
human_nomap = human[(nrow(human) - 4):nrow(human), ]
human_nomap_total = colSums(human_nomap[,2:4])
human_libsize = human_total - human_nomap_total
human = human[-c((nrow(human) - 4):nrow(human)), ]
rownames(human) = human[,1]

human0_8_exped = human[,2] > 10 | human[,3] > 10
human0_8_count = human[human0_8_exped,2:3]
colnames(human0_8_count) = c("0h","8h")
group <- factor(c(1,2))
human_deg0_8 <- DGEList(counts=human0_8_count,group=group)
# don't need to do this, because htseq-count use only uniq counted reads
#--------------------------------------------------
# y$samples$lib.size = human_libsize[1:2]
#-------------------------------------------------- 
human_deg0_8 <- calcNormFactors(human_deg0_8)
design <- model.matrix(~group)
human_deg0_8 = estimateGLMCommonDisp(human_deg0_8, method="deviance", robust=TRUE, subset=NULL)
human_fit0_8 <- glmFit(human_deg0_8, design)
human_lrt0_8 <- glmLRT(human_fit0_8,coef=2)
topTags(human_lrt0_8)
res = topTags(human_lrt0_8, n=nrow(human_lrt0_8))$table
res = cbind(rownames(res),res)
res_final = merge(res, human[,1:3],all.x=T, by.x=1, by.y=1)
colnames(res_final)[c(1,7,8)] = c("ID","Count 0h", "Count 8h")
write.table(format(res_final, digits=4),"edgeR.diff.human0_8.xls",quote=F,sep="\t",row.names=F,col.names=T)


human0_24_exped = human[,2] > 10 | human[,4] > 10
human0_24_count = human[human0_24_exped,c(2,4)]
colnames(human0_24_count) = c("0h","24h")
group <- factor(c(1,2))
human_deg0_24 <- DGEList(counts=human0_24_count,group=group)
# don't need to do this, because htseq-count use only uniq counted reads
#--------------------------------------------------
# y$samples$lib.size = human_libsize[1:2]
#-------------------------------------------------- 
human_deg0_24 <- calcNormFactors(human_deg0_24)
design <- model.matrix(~group)
human_deg0_24 = estimateGLMCommonDisp(human_deg0_24, method="deviance", robust=TRUE, subset=NULL)
human_fit0_24 <- glmFit(human_deg0_24, design)
human_lrt0_24 <- glmLRT(human_fit0_24,coef=2)
topTags(human_lrt0_24)
res = topTags(human_lrt0_24, n=nrow(human_lrt0_24))$table
res = cbind(rownames(res),res)
res_final = merge(res, human[,c(1:2,4)],all.x=T, by.x=1, by.y=1)
colnames(res_final)[c(1,7,8)] = c("ID","Count 0h", "Count 24h")
write.table(format(res_final, digits=4),"edgeR.diff.human0_24.xls",quote=F,sep="\t",row.names=F,col.names=T)

#--------------------------------------------------
# bcv <- 0.2
# counts <- matrix( rnbinom(40,size=1/bcv^2,mu=10), 20,2)
# y <- DGEList(counts=counts, group=1:2)
# et <- exactTest(y, dispersion=bcv^2)
#-------------------------------------------------- 
