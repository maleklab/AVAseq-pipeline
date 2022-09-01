library(edgeR)

args = commandArgs(trailingOnly=TRUE) 
setwd(args[1])

fileNames <- Sys.glob("*.counts")

for (fileName in fileNames) 
{
avaORF<-read.delim(fileName,header=TRUE,check.names=FALSE,stringsAsFactors = FALSE)
y<-DGEList(counts=avaORF[,2:10], genes=avaORF[,1])
ATconc<-factor(c(0,0,0,2,2,2,5,5,5))
design<-model.matrix(~ATconc)

keep <-rowSums(cpm(y)) >= 10
y<- y[keep,,keep.lib.size=FALSE]

z<-calcNormFactors(y)
rownames(design) <-colnames(z)
z <-estimateDisp(z,design,robust=TRUE)
fit <-glmFit(z,design)
lrt <-glmLRT(fit)
lrt2mM<-glmLRT(fit,coef="ATconc2")

newFileName2<-paste(sub("counts","",fileName), "2mMATdiff.txt", sep = "")
newFileName5<-paste(sub("counts","",fileName), "5mMATdiff.txt", sep = "")

write.table(topTags(lrt2mM, n=Inf),file=newFileName2,sep="\t",)
write.table(topTags(lrt, n=Inf),file=newFileName5,sep="\t",)
}
