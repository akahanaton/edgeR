library(GenomicFeatures)
txdb <- makeTxDbFromGFF("./homSap/Homo_sapiens.GRCh38.84.gff3",format="gff3")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
#--------------------------------------------------
# exonic.gene.sizes <- unlist(lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))}))
#--------------------------------------------------
exonic.gene.sizes <- unlist(sum(width(reduce(exons.list.per.gene))))

#--------------------------------------------------
# count.files = list.files("./bowtie2/", pattern="^[^J]*h.*.count2")
#--------------------------------------------------
# requir the dot (.) for th wild match
count.files = list.files("./bowtie2/", pattern="^Y.*[WC].*.count2")
#--------------------------------------------------
# count.files = count.files[-grep("JW_", count.files)]
#--------------------------------------------------

lib.size.total = read.table("./bowtie2.human.Aaron20170509.lib.size", header=F,sep=" ", stringsAsFactors=F)
lib.size.map = lib.size.total
rownames(lib.size.map) = lib.size.map[,1]
rownames(lib.size.total) = lib.size.total[,1]

i = 1
ids = vector()
for(i in 1:length(count.files)){
    id = unlist(strsplit(count.files[i],".",fixed=T))[1]
    ids[i]= id
    assign(id, read.table(paste("./bowtie2/", count.files[i],sep=""), header=F,sep="\t", stringsAsFactors=F))
    nomap.index = grep("^_", get(id)[,1])
    nomap.reads = sum(get(id)[nomap.index,2])
    lib.size.map[id,2] = lib.size.total[id,2] - nomap.reads
    assign(id, get(id)[-nomap.index,])
    assign(id, cbind(get(id), sapply(1:nrow(get(id)), function(x) {return( get(id)[x,2]/( (lib.size.map[id,2]/1000000) * (exonic.gene.sizes[get(id)[x,1]]/1000)))})))
}

datalist = lapply(1:length(ids), function(x) {return(get(ids[x]))})
rpkm.all = Reduce(function(x,y) {merge(x,y, by.x="V1", by.y="V1")}, datalist)
colnames(rpkm.all) = c("Gene",paste(c("count","rpkm"), rep(ids,each=2), sep=" "))

write.table(rpkm.all, "./rpkm.human.Aaron20170509.all.xls", quote=F, sep="\t", row.names=F, col.names=TRUE)

