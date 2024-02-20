#output_SNP_analisi

setwd("C:/Documents/outputSNPeff")
#imposta la directory nella cartella SNP
library("Biostrings")
library("seqinr")
library(Rsubread)
library(data.table)

myoutputSNP<-list.files(path=".", pattern=".VCF", all.files=FALSE, full.names=FALSE)
myoutputSNP
tabmutazioni<-data.frame(filename=myoutputSNP,sinonime=NA,upstream=NA,downstream=NA,missenso=NA)

for(aaa in 1:length(myoutputSNP))
{
a<-myoutputSNP[aaa]
#a considera scorre tra tutti i files e costituisce un nome univoco per tutti
print(a)
SNPoutput<-fread(a,skip=5)
#carica tutti i files come tabella che chiama SNPoutput
SNPoutput$Mutazioni<-NA
#crea una colonna vuota chimata mutazioni
SNPoutput$Mutazioni<-unlist(lapply(strsplit(SNPoutput$INFO,"\\|"),"[",2))
#inserisco gli elementi che vengono 2 volte dopo il| della colonna INFO nella colonna Mutazioni
sinonime<-SNPoutput$Mutazioni=="synonymous_variant"
sinonime<-sum(sinonime=="TRUE")
sinonime
tabmutazioni$sinonime[aaa]=sinonime
#mantiene solo le mutazioni sinonime e le conta
upstream<-SNPoutput$Mutazioni=="upstream_gene_variant"
upstream<-sum(upstream=="TRUE")
upstream
tabmutazioni$upstream[aaa]=upstream
#mantiene solo le mutazioni upstream e le conta
missenso<-SNPoutput$Mutazioni=="missense_variant"
missenso<-sum(missenso=="TRUE")
missenso
tabmutazioni$missenso[aaa]=missenso
#mantiene solo le mutazioni missenso e le conta
downstream<-SNPoutput$Mutazioni=="downstream_gene_variant"
downstream<-sum(downstream=="TRUE")
downstream
tabmutazioni$downstream[aaa]=downstream

#mantiene solo le downstream
}
tabmutazioni

write.table(tabmutazioni,file="tabmutazioni.txt")
#salvo una tabella in cui conta il numero per ciascuna mutazione e per ciascun file









