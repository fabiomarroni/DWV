#INPUT SNP EFF- prima esegui Freqcinquanta...
setwd("C:/Documents/Freqcinquanta")
library(data.table)
myfiles<-list.files(path=".", pattern=".VCF", all.files=FALSE, full.names=FALSE)
for(jjj in 1:length(myfiles))
{
j<-myfiles[jjj]
print(j)
SNPeff<-fread(j)
SNPeff$ALT<-SNPeff$newRef
SNPeff[ , c('cov1', 'cov2','cov3','DP','Freq.','SNP1','SNP2','SNP3','Freq2','Freq3','newRef')] <- list(NULL)
write.table(SNPeff,paste0("../inputSNPeff/",j),quote=F,row.names=F,sep="\t")
}
#serve a dare in input file SNP al programma web SNP eff insieme ad un file GTF