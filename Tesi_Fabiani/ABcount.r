#AB-APE-CONTA LA % DI READS DI VIRUS A E B MAGGIORE AL 10 %

setwd("C:/Documents/ABcount")
library(data.table)
ABcount<-fread("Virus_A_B_count.txt")
ABcount$"DWV-C"<-NULL
write.table(ABcount,"ABcount.txt",row.names=FALSE)
ABcount
#file con reads totali, reads solo di A e reads solo di B

ABcountvirusA<-ABcount[ABcount$DWVgp1>5000]
ABcountvirusA
write.table(ABcountvirusA,"ABcountvirusA.txt",row.names=FALSE)
#serve a vedere quali dei campioni di virus di APE contengono un numero di reads>5000 

ABcountvirusB<-ABcount[ABcount$VDV1_gp1>5000]
ABcountvirusB
write.table(ABcountvirusB,"ABcountvirusB.txt",row.names=FALSE)
#serve a vedere quali dei campioni di virus di APE contengono un numero di reads>5000

#SARA' DA AGGIUNERE AL VIRUS B E AL VIRUS A ALL'INIZIO



#DA ESEGUIRE SOLO 1 VOLTA
setwd("C:/Documents/etichette_virus_nuove")
getwd()
name<-fread(file="etichette 22_06.txt")
name$"NOMI SAMPLE DWV-A"<-gsub("DWV-A_","",name$"NOMI SAMPLE DWV-A")
name$"NOMI SAMPLE DWV-A"<-gsub("_viral.VCF","",name$"NOMI SAMPLE DWV-A")
ABcount$"Sample"<-gsub("ReadsPerGene.out.tab","",ABcount$"Sample")

ABcount$"Sample"[!is.na(match(ABcount$"Sample",name$"NOMI SAMPLE DWV-A"))]= na.omit(name$"NUOVE ETICHETTE DWV-A"[match(ABcount$"Sample",name$"NOMI SAMPLE DWV-A")])

ABcount$"Sample"<-gsub("DWV-A_","",ABcount$"Sample")
colnames(ABcount)<-c('Campioni','Nreads','NreadsA','NreadsB')

ABcount$"Campioni"<-gsub("","",ABcount$"Campioni")
write.table(ABcount,file="C:/Documents/ABcount/ABcount_etichettecorrette.txt", quote=FALSE,row.names=FALSE, sep="\t")
#modifica la tabella assegnando i nomi corretti ai campioni e alle colonne