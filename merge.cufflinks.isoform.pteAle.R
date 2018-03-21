 #--------------------------------------------------
 # dirs = list.files("cufflinks","^2016*")
 #--------------------------------------------------
curDir="Aaron20171208"
fq = list.files(curDir,"*.R1.*")
dirs = gsub("_R1_001.fastq.gz","", fq)

 fpkms = sapply(dirs, function(x) {
     allRpkms = read.table(paste("./cufflinks/", x, "/isoforms.fpkm_tracking",sep=""), header=T,sep="\t", stringsAsFactors=F)[,c("tracking_id","gene_id","gene_short_name","FPKM")]
     return( allRpkms )
 })

fpkms.df = matrix(unlist(fpkms), nrow=length(fpkms[[1]]),ncol=length(fpkms),byrow=F)

fpkms.final = fpkms.df[,c(1:3,1:length(dirs)*4)]
colnames(fpkms.final) = c("isoforms_id","gene_id","gene_name",dirs)

write.table(fpkms.final,paste("./rpkm.pteAle.",curDir,".cufflinks.isoforms.xls",sep=""), quote=F, sep="\t", row.names=F, col.names=TRUE)
