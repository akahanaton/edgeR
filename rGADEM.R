### R code from vignette source 'rGADEM.Rnw'
#--------------------------------------------------
# .libPaths(c(.libPaths(),"/Library/Frameworks/R.framework/Versions/3.2/Resources/library/", "/home/gmswenm/R/x86_64-pc-linux-gnu-library/3.2/"))
#--------------------------------------------------

###################################################
### code chunk number 1: loading rGADEM package
###################################################
library(rGADEM)


###################################################
### code chunk number 2: loading BSgenome package
###################################################
library(BSgenome.Hsapiens.UCSC.hg19)


###################################################
### code chunk number 3: BED File
###################################################
pwd<-"" #INPUT FILES- BedFiles, FASTA, etc.
path<- system.file("extdata/Test_100.bed",package="rGADEM")
BedFile<-paste(pwd,path,sep="")
BED<-read.table(BedFile,header=FALSE,sep="\t")
BED<-data.frame(chr=as.factor(BED[,1]),start=as.numeric(BED[,2]),end=as.numeric(BED[,3]))


###################################################
### code chunk number 4: Create the RD Files
###################################################
rgBED<-IRanges(start=BED[,2],end=BED[,3])
Sequences<-RangedData(rgBED,space=BED[,1])


###################################################
### code chunk number 5: Create the RD Files (eval = FALSE)
###################################################
##
## pwd<-"" #INPUT FILES- BedFiles, FASTA, etc.
## path<- system.file("extdata/Test_100.fasta",package="rGADEM")
## FastaFile<-paste(pwd,path,sep="")
## Sequences <- read.DNAStringSet(FastaFile, "fasta")


###################################################
### code chunk number 6: rGADEM analysis
###################################################
gadem<-GADEM(Sequences,verbose=1,genome=Hsapiens)


###################################################
### code chunk number 7: prepare PWM (eval = FALSE)
###################################################
## path<- system.file("extdata/jaspar2009.txt",package="rGADEM")
## seededPwm<-readPWMfile(path)
## grep("STAT1",names(seededPwm))
## STAT1.PWM=seededPwm[103]


###################################################
### code chunk number 8: rGADEM seeded analysis (eval = FALSE)
###################################################
## gadem<-GADEM(Sequences,verbose=1,genome=Hsapiens,Spwm=STAT1.PWM, fixSeeded=TRUE)


###################################################
### code chunk number 9: rGADEM seeded analysis (eval = FALSE)
###################################################
## gadem<-GADEM(Sequences,verbose=1,genome=Hsapiens,Spwm=STAT1.PWM)


###################################################
### code chunk number 10: pwm
###################################################
nOccurrences(gadem)


###################################################
### code chunk number 11: pwm
###################################################
nOccurrences(gadem)[1]


###################################################
### code chunk number 12: consensus
###################################################
consensus(gadem)


###################################################
### code chunk number 13: consensus
###################################################
consensus(gadem)[1]


###################################################
### code chunk number 14: position
###################################################
startPos(gadem)


###################################################
### code chunk number 15: position
###################################################
endPos(gadem)


###################################################
### code chunk number 16: parameters (eval = FALSE)
###################################################
## gadem@parameters


