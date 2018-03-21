library(edgeR)
bat0 <- read.table("./bowtie2/JW_PakiT03_Mock_0h.sam.count2", header=F,sep="\t", stringsAsFactors=F )
bat8 <- read.table("./bowtie2/JW_PakiT03_HeV_8h.sam.count2", header=F,sep="\t", stringsAsFactors=F )
bat24 <- read.table("./bowtie2/JW_PakiT03_HeV_24h.sam.count2", header=F,sep="\t", stringsAsFactors=F )
bat0_total <- read.table("./James_Transcriptome_Raw_Reads/JW_PakiT03_Mock_0h_R1.fastq.readcount", header=F,sep=" ", stringsAsFactors=F)
bat8_total <- read.table("./James_Transcriptome_Raw_Reads/JW_PakiT03_HeV_8h_R1.fastq.readcount", header=F,sep=" ", stringsAsFactors=F)
bat24_total <- read.table("./James_Transcriptome_Raw_Reads/JW_PakiT03_HeV_24h_R1.fastq.readcount", header=F,sep=" ", stringsAsFactors=F)

bat_total = c(bat0_total[,3], bat8_total[,3], bat24_total[,3])
bat = cbind(bat0, bat8[,2], bat24[,2])
bat_nomap = bat[(nrow(bat) - 4):nrow(bat), ]
bat_nomap_total = colSums(bat_nomap[,2:4])
bat_libsize = bat_total - bat_nomap_total
bat = bat[-c((nrow(bat) - 4):nrow(bat)), ]
rownames(bat) = bat[,1]

bat0_8_exped = bat[,2] > 10 | bat[,3] > 10
bat0_8_count = bat[bat0_8_exped,2:3]
colnames(bat0_8_count) = c("0h","8h")
group <- factor(c(1,2))
bat_deg0_8 <- DGEList(counts=bat0_8_count,group=group)
# don't need to do this, because htseq-count use only uniq counted reads
#--------------------------------------------------
# y$samples$lib.size = bat_libsize[1:2]
#-------------------------------------------------- 
bat_deg0_8 <- calcNormFactors(bat_deg0_8)
design <- model.matrix(~group)
bat_deg0_8 = estimateGLMCommonDisp(bat_deg0_8, method="deviance", robust=TRUE, subset=NULL)
bat_fit0_8 <- glmFit(bat_deg0_8, design)
bat_lrt0_8 <- glmLRT(bat_fit0_8,coef=2)
topTags(bat_lrt0_8)
res = topTags(bat_lrt0_8, n=nrow(bat_lrt0_8))$table
res = cbind(rownames(res),res)
res_final = merge(res, bat[,1:3],all.x=T, by.x=1, by.y=1)
colnames(res_final)[c(1,7,8)] = c("ID","Count 0h", "Count 8h")
write.table(format(res_final, digits=4),"edgeR.diff.bat0_8.xls",quote=F,sep="\t",row.names=F,col.names=T)


bat0_24_exped = bat[,2] > 10 | bat[,4] > 10
bat0_24_count = bat[bat0_24_exped,c(2,4)]
colnames(bat0_24_count) = c("0h","24h")
group <- factor(c(1,2))
bat_deg0_24 <- DGEList(counts=bat0_24_count,group=group)
# don't need to do this, because htseq-count use only uniq counted reads
#--------------------------------------------------
# y$samples$lib.size = bat_libsize[1:2]
#-------------------------------------------------- 
bat_deg0_24 <- calcNormFactors(bat_deg0_24)
design <- model.matrix(~group)
bat_deg0_24 = estimateGLMCommonDisp(bat_deg0_24, method="deviance", robust=TRUE, subset=NULL)
bat_fit0_24 <- glmFit(bat_deg0_24, design)
bat_lrt0_24 <- glmLRT(bat_fit0_24,coef=2)
topTags(bat_lrt0_24)
res = topTags(bat_lrt0_24, n=nrow(bat_lrt0_24))$table
res = cbind(rownames(res),res)
res_final = merge(res, bat[,c(1:2,4)],all.x=T, by.x=1, by.y=1)
colnames(res_final)[c(1,7,8)] = c("ID","Count 0h", "Count 24h")
write.table(format(res_final, digits=4),"edgeR.diff.bat0_24.xls",quote=F,sep="\t",row.names=F,col.names=T)
