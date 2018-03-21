 dirs = list.files("cufflinks","^J.*ES*")

 id = dirs[1]
 allRpkms = read.table(paste("./cufflinks/", id, "/genes.fpkm_tracking",sep=""), header=T,sep="\t", stringsAsFactors=F)
 allRpkms = allRpkms[-grep("CUFF",allRpkms[,1]),c("gene_id","FPKM")]

 id=2
 for (id in 2:length(dirs)){
     tmp = read.table(paste("./cufflinks/", dirs[id], "/genes.fpkm_tracking",sep=""), header=T,sep="\t", stringsAsFactors=F)
     tmp = tmp[-grep("CUFF",tmp[,1]),]
     allRpkms = merge(allRpkms,tmp[,c("gene_id","FPKM")],by.x = "gene_id",by.y="gene_id",all=T)
 }
 allRpkms[is.na(allRpkms)] = 0
 #--------------------------------------------------
 # fpkm_idx = grep("FPKM",colnames(allRpkms))
 # name_idx = grep("name",colnames(allRpkms))
 # gene_name = sapply(1:nrow(allRpkms), function(x) { name = allRpkms[x,name_idx]; name = name[name != 0]; return(unique(name)) })
 # allRpkms.final = cbind(allRpkms[,1],gene_name, allRpkms[,fpkm_idx])
 #--------------------------------------------------
 colnames(allRpkms) = c("Gene",dirs)

write.table(allRpkms, "./rpkm.eonSpe.Aaron20170412.ES.tophat.xls", quote=F, sep="\t", row.names=F, col.names=TRUE)
