library(edgeR)
library(GenomicFeatures)
count.files =list.files("./htseq","^2016*")

lib.size.total = read.table("./bowtie2.Aaron20170904.lib.size", header=F,sep=" ", stringsAsFactors=F)
lib.size.map = lib.size.total
rownames(lib.size.map) = lib.size.map[,1]
rownames(lib.size.total) = lib.size.total[,1]

exp.cutoff = 10

#--------------------------------------------------
# txdb <- makeTxDbFromGFF("./pteAle/GCF_000325575.1_ASM32557v1_genomic.gff",format="gff3")
#--------------------------------------------------
#--------------------------------------------------
# txdb <- makeTxDbFromGFF("/synology/Ming/genome_assembly/P_alecto/GCF_000325575.1_ASM32557v1_genomic.gff",format="gff3")
#--------------------------------------------------
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene",use.names=F)
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
#--------------------------------------------------
# exonic.gene.sizes <- unlist(lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))}))
#--------------------------------------------------
exonic.gene.sizes <- unlist(sum(width(reduce(exons.list.per.gene))))

i = 1
ids = vector()
for(i in 1:length(count.files)){
    id = unlist(strsplit(count.files[i],".",fixed=T))[1]
    ids[i]= id
    assign(id, read.table(paste("./htseq/", count.files[i],sep=""), header=F,sep="\t", stringsAsFactors=F))
    nomap.index = grep("^_", get(id)[,1])
    nomap.reads = sum(get(id)[nomap.index,2])
    lib.size.map[id,2] = lib.size.total[id,2] - nomap.reads
    assign(id, get(id)[-nomap.index,]) # remove the read statistical lines in the reads count results
    assign(id, cbind(get(id), sapply(1:nrow(get(id)), function(x) {return( (get(id)[x,2]+1)/( (lib.size.map[id,2]/1000000) * (exonic.gene.sizes[get(id)[x,1]]/1000)))})))
}
datalist = lapply(1:length(ids), function(x) {return(get(ids[x]))})
rpkm.all = Reduce(function(x,y) {merge(x,y, by.x="V1", by.y="V1")}, datalist)
colnames(rpkm.all) = c("Gene",paste(c("","rpkm "), rep(ids,each=2), sep=""))
rownames(rpkm.all) = rpkm.all[,1]

group <- factor(c(rep(1,5),rep(2,5)))

liver.index = grep("Liver",ids,value=T)
liver.lib = lib.size.map[liver.index,]
liver.df = rpkm.all[,liver.index]
liver.df.exped = liver.df > exp.cutoff
liver.index.control = seq(1, length(liver.index), 3)
liver.index.lps = seq(1, length(liver.index), 3) + 1
liver.df.exped.lps = rowSums(liver.df.exped[,liver.index.control]) > 3 | rowSums(liver.df.exped[,liver.index.lps]) > 3
liver.df.lps = liver.df[liver.df.exped.lps, c(liver.index.control, liver.index.lps)]
liver.lib.lps = liver.lib[c(liver.index.control, liver.index.lps),]
liver.deglist.lps <- DGEList(counts=liver.df.lps,group=group)
liver.deglist.lps$samples$lib.size = liver.lib.lps[,2]
liver.deglist.lps <- calcNormFactors(liver.deglist.lps)
design <- model.matrix(~group)
liver.deglist.lps = estimateGLMCommonDisp(liver.deglist.lps, method="deviance", robust=TRUE, subset=NULL)
liver.deglist.lps.fit <- glmFit(liver.deglist.lps, design)
liver.deglist.lps.lrt = glmLRT(liver.deglist.lps.fit, coef=2)
liver.deglist.lps.res = topTags(liver.deglist.lps.lrt, n=nrow(liver.deglist.lps.lrt))$table
liver.deglist.lps.res = cbind(rownames(liver.deglist.lps.res),liver.deglist.lps.res)
liver.index.pic = seq(1, length(liver.index), 3) + 2
liver.df.exped.pic = rowSums(liver.df.exped[,liver.index.control]) > 3 | rowSums(liver.df.exped[,liver.index.pic]) > 3
liver.df.pic = liver.df[liver.df.exped.pic, c(liver.index.control, liver.index.pic)]
liver.lib.pic = liver.lib[c(liver.index.control, liver.index.pic),]
liver.deglist.pic <- DGEList(counts=liver.df.pic,group=group)
liver.deglist.pic$samples$lib.size = liver.lib.pic[,2]
liver.deglist.pic <- calcNormFactors(liver.deglist.pic)
design <- model.matrix(~group)
liver.deglist.pic = estimateGLMCommonDisp(liver.deglist.pic, method="deviance", robust=TRUE, subset=NULL)
liver.deglist.pic.fit <- glmFit(liver.deglist.pic, design)
liver.deglist.pic.lrt = glmLRT(liver.deglist.pic.fit, coef=2)
liver.deglist.pic.res = topTags(liver.deglist.pic.lrt, n=nrow(liver.deglist.pic.lrt))$table
liver.deglist.pic.res = cbind(rownames(liver.deglist.pic.res),liver.deglist.pic.res)
liver.rpkm = cbind(rpkm.all[,c("Gene",grep("Liver", colnames(rpkm.all),value=T))], liver.df.exped.lps, liver.df.exped.pic)
liver.res.list = list(liver.rpkm, liver.deglist.lps.res, liver.deglist.pic.res)
liver.res = Reduce(function(x, y) merge(x, y, by.x=1,by.y=1, all=TRUE), liver.res.list, accumulate=FALSE)
write.table(liver.res,file="Aaron20170904.liver.xls", quote=F,sep="\t",col.names=T,row.names=F)

ln.index = grep("LN",ids,value=T)
ln.lib = lib.size.map[ln.index,]
ln.df = rpkm.all[,ln.index]
ln.df.exped = ln.df > exp.cutoff
ln.index.control = seq(1, length(ln.index), 3)
ln.index.lps = seq(1, length(ln.index), 3) + 1
ln.df.exped.lps = rowSums(ln.df.exped[,ln.index.control]) > 3 | rowSums(ln.df.exped[,ln.index.lps]) > 3
ln.df.lps = ln.df[ln.df.exped.lps, c(ln.index.control, ln.index.lps)]
ln.lib.lps = ln.lib[c(ln.index.control, ln.index.lps),]
ln.deglist.lps <- DGEList(counts=ln.df.lps,group=group)
ln.deglist.lps$samples$lib.size = ln.lib.lps[,2]
ln.deglist.lps <- calcNormFactors(ln.deglist.lps)
design <- model.matrix(~group)
ln.deglist.lps = estimateGLMCommonDisp(ln.deglist.lps, method="deviance", robust=TRUE, subset=NULL)
ln.deglist.lps.fit <- glmFit(ln.deglist.lps, design)
ln.deglist.lps.lrt = glmLRT(ln.deglist.lps.fit, coef=2)
ln.deglist.lps.res = topTags(ln.deglist.lps.lrt, n=nrow(ln.deglist.lps.lrt))$table
ln.deglist.lps.res = cbind(rownames(ln.deglist.lps.res),ln.deglist.lps.res)
ln.index.pic = seq(1, length(ln.index), 3) + 2
ln.df.exped.pic = rowSums(ln.df.exped[,ln.index.control]) > 3 | rowSums(ln.df.exped[,ln.index.pic]) > 3
ln.df.pic = ln.df[ln.df.exped.pic, c(ln.index.control, ln.index.pic)]
ln.lib.pic = ln.lib[c(ln.index.control, ln.index.pic),]
ln.deglist.pic <- DGEList(counts=ln.df.pic,group=group)
ln.deglist.pic$samples$lib.size = ln.lib.pic[,2]
ln.deglist.pic <- calcNormFactors(ln.deglist.pic)
design <- model.matrix(~group)
ln.deglist.pic = estimateGLMCommonDisp(ln.deglist.pic, method="deviance", robust=TRUE, subset=NULL)
ln.deglist.pic.fit <- glmFit(ln.deglist.pic, design)
ln.deglist.pic.lrt = glmLRT(ln.deglist.pic.fit, coef=2)
ln.deglist.pic.res = topTags(ln.deglist.pic.lrt, n=nrow(ln.deglist.pic.lrt))$table
ln.deglist.pic.res = cbind(rownames(ln.deglist.pic.res),ln.deglist.pic.res)
ln.rpkm = cbind(rpkm.all[,c("Gene",grep("LN", colnames(rpkm.all),value=T))], ln.df.exped.lps, ln.df.exped.pic)
ln.res.list = list(ln.rpkm, ln.deglist.lps.res, ln.deglist.pic.res)
ln.res = Reduce(function(x, y) merge(x, y, by.x=1,by.y=1, all=TRUE), ln.res.list, accumulate=FALSE)
write.table(ln.res,file="Aaron20170904.ln.xls", quote=F,sep="\t",col.names=T,row.names=F)

lung.index = grep("Lung",ids,value=T)
lung.lib = lib.size.map[lung.index,]
lung.df = rpkm.all[,lung.index]
lung.df.exped = lung.df > exp.cutoff
lung.index.control = seq(1, length(lung.index), 3)
lung.index.lps = seq(1, length(lung.index), 3) + 1
lung.df.exped.lps = rowSums(lung.df.exped[,lung.index.control]) > 3 | rowSums(lung.df.exped[,lung.index.lps]) > 3
lung.df.lps = lung.df[lung.df.exped.lps, c(lung.index.control, lung.index.lps)]
lung.lib.lps = lung.lib[c(lung.index.control, lung.index.lps),]
lung.deglist.lps <- DGEList(counts=lung.df.lps,group=group)
lung.deglist.lps$samples$lib.size = lung.lib.lps[,2]
lung.deglist.lps <- calcNormFactors(lung.deglist.lps)
design <- model.matrix(~group)
lung.deglist.lps = estimateGLMCommonDisp(lung.deglist.lps, method="deviance", robust=TRUE, subset=NULL)
lung.deglist.lps.fit <- glmFit(lung.deglist.lps, design)
lung.deglist.lps.lrt = glmLRT(lung.deglist.lps.fit, coef=2)
lung.deglist.lps.res = topTags(lung.deglist.lps.lrt, n=nrow(lung.deglist.lps.lrt))$table
lung.deglist.lps.res = cbind(rownames(lung.deglist.lps.res),lung.deglist.lps.res)
lung.index.pic = seq(1, length(lung.index), 3) + 2
lung.df.exped.pic = rowSums(lung.df.exped[,lung.index.control]) > 3 | rowSums(lung.df.exped[,lung.index.pic]) > 3
lung.df.pic = lung.df[lung.df.exped.pic, c(lung.index.control, lung.index.pic)]
lung.lib.pic = lung.lib[c(lung.index.control, lung.index.pic),]
lung.deglist.pic <- DGEList(counts=lung.df.pic,group=group)
lung.deglist.pic$samples$lib.size = lung.lib.pic[,2]
lung.deglist.pic <- calcNormFactors(lung.deglist.pic)
design <- model.matrix(~group)
lung.deglist.pic = estimateGLMCommonDisp(lung.deglist.pic, method="deviance", robust=TRUE, subset=NULL)
lung.deglist.pic.fit <- glmFit(lung.deglist.pic, design)
lung.deglist.pic.lrt = glmLRT(lung.deglist.pic.fit, coef=2)
lung.deglist.pic.res = topTags(lung.deglist.pic.lrt, n=nrow(lung.deglist.pic.lrt))$table
lung.deglist.pic.res = cbind(rownames(lung.deglist.pic.res),lung.deglist.pic.res)
lung.rpkm = cbind(rpkm.all[,c("Gene",grep("Lung", colnames(rpkm.all),value=T))], lung.df.exped.lps, lung.df.exped.pic)
lung.res.list = list(lung.rpkm, lung.deglist.lps.res, lung.deglist.pic.res)
lung.res = Reduce(function(x, y) merge(x, y, by.x=1,by.y=1, all=TRUE), lung.res.list, accumulate=FALSE)
write.table(lung.res,file="Aaron20170904.lung.xls", quote=F,sep="\t",col.names=T,row.names=F)

spleen.index = grep("Spleen",ids,value=T)
spleen.lib = lib.size.map[spleen.index,]
spleen.df = rpkm.all[,spleen.index]
spleen.df.exped = spleen.df > exp.cutoff
spleen.index.control = seq(1, length(spleen.index), 3)
spleen.index.lps = seq(1, length(spleen.index), 3) + 1
spleen.df.exped.lps = rowSums(spleen.df.exped[,spleen.index.control]) > 3 | rowSums(spleen.df.exped[,spleen.index.lps]) > 3
spleen.df.lps = spleen.df[spleen.df.exped.lps, c(spleen.index.control, spleen.index.lps)]
spleen.lib.lps = spleen.lib[c(spleen.index.control, spleen.index.lps),]
spleen.deglist.lps <- DGEList(counts=spleen.df.lps,group=group)
spleen.deglist.lps$samples$lib.size = spleen.lib.lps[,2]
spleen.deglist.lps <- calcNormFactors(spleen.deglist.lps)
design <- model.matrix(~group)
spleen.deglist.lps = estimateGLMCommonDisp(spleen.deglist.lps, method="deviance", robust=TRUE, subset=NULL)
spleen.deglist.lps.fit <- glmFit(spleen.deglist.lps, design)
spleen.deglist.lps.lrt = glmLRT(spleen.deglist.lps.fit, coef=2)
spleen.deglist.lps.res = topTags(spleen.deglist.lps.lrt, n=nrow(spleen.deglist.lps.lrt))$table
spleen.deglist.lps.res = cbind(rownames(spleen.deglist.lps.res),spleen.deglist.lps.res)
spleen.index.pic = seq(1, length(spleen.index), 3) + 2
spleen.df.exped.pic = rowSums(spleen.df.exped[,spleen.index.control]) > 3 | rowSums(spleen.df.exped[,spleen.index.pic]) > 3
spleen.df.pic = spleen.df[spleen.df.exped.pic, c(spleen.index.control, spleen.index.pic)]
spleen.lib.pic = spleen.lib[c(spleen.index.control, spleen.index.pic),]
spleen.deglist.pic <- DGEList(counts=spleen.df.pic,group=group)
spleen.deglist.pic$samples$lib.size = spleen.lib.pic[,2]
spleen.deglist.pic <- calcNormFactors(spleen.deglist.pic)
design <- model.matrix(~group)
spleen.deglist.pic = estimateGLMCommonDisp(spleen.deglist.pic, method="deviance", robust=TRUE, subset=NULL)
spleen.deglist.pic.fit <- glmFit(spleen.deglist.pic, design)
spleen.deglist.pic.lrt = glmLRT(spleen.deglist.pic.fit, coef=2)
spleen.deglist.pic.res = topTags(spleen.deglist.pic.lrt, n=nrow(spleen.deglist.pic.lrt))$table
spleen.deglist.pic.res = cbind(rownames(spleen.deglist.pic.res),spleen.deglist.pic.res)
spleen.rpkm = cbind(rpkm.all[,c("Gene",grep("Spleen", colnames(rpkm.all),value=T))], spleen.df.exped.lps, spleen.df.exped.pic)
spleen.res.list = list(spleen.rpkm, spleen.deglist.lps.res, spleen.deglist.pic.res)
spleen.res = Reduce(function(x, y) merge(x, y, by.x=1,by.y=1, all=TRUE), spleen.res.list, accumulate=FALSE)
write.table(spleen.res,file="Aaron20170904.spleen.xls", quote=F,sep="\t",col.names=T,row.names=F)
