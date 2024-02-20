#STATISTICHE PER GENE:per ogni campione e per ogni gene individua il numero di mutazioni per n° di basi azotate - 17/05

setwd("C:/Documents/outputSNPeff")
#imposta la directory nella cartella outputSNPeff
library("Biostrings")
library("seqinr")
library(Rsubread)
library(data.table)

sinomiss<- data.frame(
GENE = c("Lp","VP3","VP4","VP1","VP2","Elicasi","VPg","3C+RdRp"),
POSi = c("1145","1779","2518","2582","3831","4994","7411","7673"),
POSf = c("1778","2517","2581","3830","4605","6410","7672","9812"),stringsAsFactors = FALSE)
#crea un dataframe in cui ci sono i nomi dei geni e le loro posizioni sulla base del "DALMON"

myfiles<-list.files(pattern=".VCF")
for (aaa in 1:length (myfiles))
{
elab<-fread(myfiles[aaa],skip=5)
#leggo un file qualsiasi degli outputsnp
elab$Mutazioni<-NA
#crea una colonna vuota chimata mutazioni
elab$Mutazioni<-unlist(lapply(strsplit(elab$INFO,"\\|"),"[",2))
#prendo il secondo valore dopo | e lo incollo nella nuova colonna mutazioni

Lpsin<-elab[elab$POS>=1145 & elab$POS<=1778 & elab$Mutazioni=="synonymous_variant",]
#sulla base delle posizioni identificate dal gene LP cerco le sole mutazioni sinonime
sinomiss$nSINONIME<-NA
sinomiss$nSINONIME[1]<-nrow(Lpsin)
#inserisco il numero nel dataframe di partenza

Lpmiss<-elab[elab$POS>=1145 & elab$POS<=1778 & elab$Mutazioni=="missense_variant",]
#sulla base delle posizioni identificate dal gene LP cerco le sole mutazioni missenso
sinomiss$nMISSENSO<-NA
sinomiss$nMISSENSO[1]<-nrow(Lpmiss)


VP3sin<-elab[elab$POS>=1779 & elab$POS<=2517 & elab$Mutazioni=="synonymous_variant",]
#sulla base delle posizioni identificate dal gene VP3 cerco le sole mutazioni sinonime
sinomiss$nSINONIME[2]<-nrow(VP3sin)

VP3miss<-elab[elab$POS>=1779 & elab$POS<=2517 & elab$Mutazioni=="missense_variant",]
#sulla base delle posizioni identificate dal gene VP3 cerco le sole mutazioni missenso
sinomiss$nMISSENSO[2]<-nrow(VP3miss)


VP4sin<-elab[elab$POS>=2518 & elab$POS<=2581 & elab$Mutazioni=="synonymous_variant",]
#sulla base delle posizioni identificate dal gene VP4 cerco le sole mutazioni sinonime
sinomiss$nSINONIME[3]<-nrow(VP4sin)

VP4miss<-elab[elab$POS>=2518 & elab$POS<=2581 & elab$Mutazioni=="missense_variant",]
#sulla base delle posizioni identificate dal gene VP4 cerco le sole mutazioni missenso
sinomiss$nMISSENSO[3]<-nrow(VP4miss)


VP1sin<-elab[elab$POS>=2582 & elab$POS<=3830 & elab$Mutazioni=="synonymous_variant",]
#sulla base delle posizioni identificate dal gene VP1 cerco le sole mutazioni sinonime
sinomiss$nSINONIME[4]<-nrow(VP1sin)

VP1miss<-elab[elab$POS>=2582 & elab$POS<=3830 & elab$Mutazioni=="missense_variant",]
#sulla base delle posizioni identificate dal gene VP1 cerco le sole mutazioni missenso
sinomiss$nMISSENSO[4]<-nrow(VP1miss)


VP2sin<-elab[elab$POS>=3831 & elab$POS<=4605 & elab$Mutazioni=="synonymous_variant",]
#sulla base delle posizioni identificate dal gene VP2 cerco le sole mutazioni sinonime
sinomiss$nSINONIME[5]<-nrow(VP2sin)

VP2miss<-elab[elab$POS>=3831 & elab$POS<=4605 & elab$Mutazioni=="missense_variant",]
#sulla base delle posizioni identificate dal gene VP2 cerco le sole mutazioni missenso
sinomiss$nMISSENSO[5]<-nrow(VP2miss)


elicasisin<-elab[elab$POS>=4994 & elab$POS<=6410 & elab$Mutazioni=="synonymous_variant",]
#sulla base delle posizioni identificate dal gene elicasi cerco le sole mutazioni sinonime
sinomiss$nSINONIME[6]<-nrow(elicasisin)

elicasimiss<-elab[elab$POS>=4994 & elab$POS<=6410 & elab$Mutazioni=="missense_variant",]
#sulla base delle posizioni identificate dal gene elicasi cerco le sole mutazioni missenso
sinomiss$nMISSENSO[6]<-nrow(elicasimiss)


VPgsin<-elab[elab$POS>=7411 & elab$POS<=7672 & elab$Mutazioni=="synonymous_variant",]
#sulla base delle posizioni identificate dal gene VPg cerco le sole mutazioni sinonime
sinomiss$nSINONIME[7]<-nrow(VPgsin)

VPgmiss<-elab[elab$POS>=7411 & elab$POS<=7672 & elab$Mutazioni=="missense_variant",]
#sulla base delle posizioni identificate dal gene VPg cerco le sole mutazioni missenso
sinomiss$nMISSENSO[7]<-nrow(VPgmiss)


treCRdRpsin<-elab[elab$POS>=7673 & elab$POS<=9812 & elab$Mutazioni=="synonymous_variant",]
#sulla base delle posizioni identificate dal gene 3C e RdRp cerco le sole mutazioni sinonime
sinomiss$nSINONIME[8]<-nrow(treCRdRpsin)

treCRdRpmiss<-elab[elab$POS>=7673 & elab$POS<=9812 & elab$Mutazioni=="missense_variant",]
#sulla base delle posizioni identificate dal gene 3C e RdRp cerco le sole mutazioni missenso
sinomiss$nMISSENSO[8]<-nrow(treCRdRpmiss)

sum(sinomiss$nSINONIME)
#somma tutte le mutazioni sinonime dei vari geni di un singolo campione di virus
sum(sinomiss$nMISSENSO)
#somma tutte le mutazioni missenso dei vari geni di un singolo campione di virus

sinomiss$GenedivS<-NA
sinomiss$GenedivS[1]<-sinomiss$nSINONIME[1]/(1778-1145)
#frequenza sinonime per gene LP
sinomiss$GenedivS[2]<-sinomiss$nSINONIME[2]/(2517-1779)
#frequenza sinonime per gene VP3
sinomiss$GenedivS[3]<-sinomiss$nSINONIME[3]/(2581-2518)
#frequenza sinonime per gene VP4
sinomiss$GenedivS[4]<-sinomiss$nSINONIME[4]/(3830-2582)
#frequenza sinonime per gene VP1
sinomiss$GenedivS[5]<-sinomiss$nSINONIME[5]/(4605-3831)
#frequenza sinonime per gene VP2
sinomiss$GenedivS[6]<-sinomiss$nSINONIME[6]/(6410-4994)
#frequenza sinonime per gene elicasi
sinomiss$GenedivS[7]<-sinomiss$nSINONIME[7]/(7672-7411)
#frequenza sinonime per gene VPg
sinomiss$GenedivS[8]<-sinomiss$nSINONIME[8]/(9812-7673)
#frequenza sinonime per gene 3C+RdRp
sinomiss$GenedivS<-round(sinomiss$GenedivS,digit=4)
#arrotondo la freq al 4 numero decimale
sinomiss



sinomiss$GenedivM<-NA
sinomiss$GenedivM[1]<-sinomiss$nMISSENSO[1]/(1778-1145)
#frequenza missenso per gene LP
sinomiss$GenedivM[2]<-sinomiss$nMISSENSO[2]/(2517-1779)
#frequenza missenso per gene VP3
sinomiss$GenedivM[3]<-sinomiss$nMISSENSO[3]/(2581-2518)
#frequenza missenso per gene VP4
sinomiss$GenedivM[4]<-sinomiss$nMISSENSO[4]/(3830-2582)
#frequenza missenso per gene VP1
sinomiss$GenedivM[5]<-sinomiss$nMISSENSO[5]/(4605-3831)
#frequenza missenso per gene VP2
sinomiss$GenedivM[6]<-sinomiss$nMISSENSO[6]/(6410-4994)
#frequenza missenso per gene elicasi
sinomiss$GenedivM[7]<-sinomiss$nMISSENSO[7]/(7672-7411)
#frequenza missenso per gene VPg
sinomiss$GenedivM[8]<-sinomiss$nMISSENSO[8]/(9812-7673)
#frequenza missenso per gene 3C+RdRp
sinomiss$GenedivM<-round(sinomiss$GenedivM,digit=4)
#arrotondo la freq al 4 numero decimale
sinomiss

outfile<-gsub(".VCF",".txt",myfiles[aaa])
outfile<-paste0("sinomiss_",outfile)
outfile<-paste0("C:/Documents/mutazioni_Sinmiss/",outfile)
write.table(sinomiss, file=outfile, quote=FALSE,row.names=FALSE,sep="\t")

}
#mi salva un file di testo per ciascun campione di A e di B in cui si riportano n° di mutazioni missenso e consenso per ciascun gene e anche la frequenza per gene (al momento non in %)


#TABELLA TOT SINONIME, TOT MISSENSO, TOTgenedivpoliproteinaS e M PER CIASCUN CAMPIONE DI DWV-A - 15/06

setwd("C:/Documents/mutazioni_Sinmiss")
#imposta la directory nella cartella mutazioni Sinmiss

letturacampioniA<-list.files(pattern="DWV-A")
letturacampioniA
#lista solo i virus A

totA<- data.frame(Campione = letturacampioniA, Sinonimetot=NA, Missensotot=NA, GenedivStot=NA, GenedivMtot=NA,stringsAsFactors=FALSE)

for (aaa in 1:length (letturacampioniA))
{
a<-letturacampioniA[aaa]
leggitabella<-fread(a,data.table=FALSE)


totA$Sinonimetot[aaa]<-sum(leggitabella$nSINONIME)
totA$Missensotot[aaa]<-sum(leggitabella$nMISSENSO)
totA$GenedivStot[aaa]<-sum(leggitabella$GenedivS)
totA$GenedivMtot[aaa]<-sum(leggitabella$GenedivM)
}

setwd("C:/Documents/etichette_virus_nuove")
getwd()
namesA<-fread(file="etichette api tesi_MF.txt")
namesA
totA$Campione<-gsub("sinomiss_","",totA$Campione)
totA$Campione<-gsub("_snpeff.txt",".VCF",totA$Campione)
totA$Campione[!is.na(match(totA$Campione,namesA$"NOMI SAMPLE DWV-A"))]= na.omit(namesA$"NUOVE ETICHETTE DWV-A"[match(totA$Campione,namesA$"NOMI SAMPLE DWV-A")])
totA$Campione<-gsub("DWV-A_","",totA$Campione)
#modificare nomi etichette

write.table(totA, file="C:/Documents/mutazioni_Sinmiss/totpoliproteinaA.txt", quote=FALSE, row.names=FALSE, sep="\t")


#TABELLA TOT SINONIME E TOT MISSENSO PER CIASCUN CAMPIONE DI DWV-B - 15/06

setwd("C:/Documents/mutazioni_Sinmiss")
#imposta la directory nella cartella mutazioni Sinmiss

letturacampioniB<-list.files(pattern="DWV-B")
letturacampioniB
#lista solo i virus A

totB<- data.frame(Campione = letturacampioniB, Sinonimetot=NA, Missensotot=NA, GenedivStot=NA, GenedivMtot=NA,stringsAsFactors=FALSE)

for (aaa in 1:length (letturacampioniB))
{
a<-letturacampioniB[aaa]
leggitabella<-fread(a,data.table=FALSE)

totB$Sinonimetot[aaa]<-sum(leggitabella$nSINONIME)
totB$Missensotot[aaa]<-sum(leggitabella$nMISSENSO)
totB$GenedivStot[aaa]<-sum(leggitabella$GenedivS)
totB$GenedivMtot[aaa]<-sum(leggitabella$GenedivM)
}

setwd("C:/Documents/etichette_virus_nuove")
getwd()
namesB<-fread(file="etichette api tesi_MF.txt")
namesB
totB$Campione<-gsub("sinomiss_","",totB$Campione)
totB$Campione<-gsub("_snpeff.txt",".VCF",totB$Campione)
totB$Campione[!is.na(match(totB$Campione,namesB$"NOMI SAMPLE DWV-B"))]= na.omit(namesB$"NUOVE ETICHETTE DWV-B"[match(totB$Campione,namesB$"NOMI SAMPLE DWV-B")])
totB$Campione<-gsub("DWV-B_","",totB$Campione)
#modifica nomi etichette

write.table(totB, file="C:/Documents/mutazioni_Sinmiss/totpoliproteinaB.txt", quote=FALSE, row.names=FALSE, sep="\t")



#MEDIA MUTAZIONI PER GENE PER DWV-A - 23/05

setwd("C:/Documents/mutazioni_Sinmiss")
#imposta la directory nella cartella mutazioni Sinmiss

mysinomiss<-list.files(pattern="DWV-A")
mysinomiss
#lista solo i virus A

mydfA<- data.frame(
SAMPLE = mysinomiss,LpS=NA,LpM=NA,VP3S=NA,VP3M=NA,VP4S=NA,VP4M=NA,VP1S=NA,VP1M=NA,VP2S=NA,VP2M=NA,
ElicasiS=NA,ElicasiM=NA,VPgS=NA,VPgM=NA,treCRdRpS=NA,treCRdRpM=NA,stringsAsFactors=FALSE)

for (aaa in 1:nrow (mydfA))
{
a<-mysinomiss[aaa]
totalsinomiss<-fread(a,data.table=FALSE)

mydfA$LpS[aaa]<-totalsinomiss[1,6]
mydfA$LpM[aaa]<-totalsinomiss[1,7]
mydfA$VP3S[aaa]<-totalsinomiss[2,6]
mydfA$VP3M[aaa]<-totalsinomiss[2,7]
mydfA$VP4S[aaa]<-totalsinomiss[3,6]
mydfA$VP4M[aaa]<-totalsinomiss[3,7]
mydfA$VP1S[aaa]<-totalsinomiss[4,6]
mydfA$VP1M[aaa]<-totalsinomiss[4,7]
mydfA$VP2S[aaa]<-totalsinomiss[5,6]
mydfA$VP2M[aaa]<-totalsinomiss[5,7]
mydfA$ElicasiS[aaa]<-totalsinomiss[6,6]
mydfA$ElicasiM[aaa]<-totalsinomiss[6,7]
mydfA$VPgS[aaa]<-totalsinomiss[7,6]
mydfA$VPgM[aaa]<-totalsinomiss[7,7]
mydfA$treCRdRpS[aaa]<-totalsinomiss[8,6]
mydfA$treCRdRpM[aaa]<-totalsinomiss[8,7]
}

write.table(mydfA, file="C:/Documents/mutazioni_Sinmiss/allgenesDWVA.txt", quote=FALSE, row.names=FALSE, sep="\t")



#GRAFICI E TEST VIRUS A SINONIME PER TUTTI I GENI - 03/06

setwd("C:/Documents/png")
png("boxplot_sinonime_DWV-A.png",units="cm",res=600, width=15, height=15)
virus_A_sin_allgenes<-boxplot(mydfA$LpS,mydfA$VP3S,mydfA$VP4S,mydfA$VP1S,mydfA$VP2S,mydfA$ElicasiS,mydfA$VPgS,mydfA$treCRdRpS, main="Boxplot mutazioni sinonime DWV-A", names=c("LpS","VP3S","VP4S","VP1S","VP2S","ElicasiS","VPgS","3C-RdRpS"), ylab="Gene diversity",las=2, ylim=c(0,0.05)) 
virus_A_sin_allgenes
dev.off()
#salvo le mutazioni sinonime dei virus A in un file png


#GRAFICI E TEST VIRUS A MISSENSO PER TUTTI I GENI - 03/06

png("boxplot_missenso_DWV-A.png",units="cm",res=600, width=15, height=15)
virus_A_miss_allgenes<-boxplot(mydfA$LpM,mydfA$VP3M,mydfA$VP4M,mydfA$VP1M,mydfA$VP2M,mydfA$ElicasiM,mydfA$VPgM,mydfA$treCRdRpM, main="Boxplot mutazioni missenso DWV-A", names=c("LpM","VP3M","VP4M","VP1M","VP2M","ElicasiM","VPgM","3C-RdRpM"), ylab="Gene diversity",las=2, ylim=c(0,0.03))
virus_A_miss_allgenes
dev.off()



#MEDIA MUTAZIONI PER GENE PER DWV-B - 23/05

setwd("C:/Documents/mutazioni_Sinmiss")
#imposta la directory nella cartella mutazioni Sinmiss

mysinomiss<-list.files(pattern="DWV-B")
mysinomiss
#lista solo i virus B

mydfB<- data.frame(
SAMPLE = mysinomiss,LpS=NA,LpM=NA,VP3S=NA,VP3M=NA,VP4S=NA,VP4M=NA,VP1S=NA,VP1M=NA,VP2S=NA,VP2M=NA,
ElicasiS=NA,ElicasiM=NA,VPgS=NA,VPgM=NA,treCRdRpS=NA,treCRdRpM=NA,stringsAsFactors=FALSE)


for (aaa in 1:nrow (mydfB))
{
a<-mysinomiss[aaa]
totalsinomiss<-fread(a,data.table=FALSE)
#totalsinoLp<-totalsinomiss[totalsinomiss$FreqS[1],]
mydfB$LpS[aaa]<-totalsinomiss[1,6]
mydfB$LpM[aaa]<-totalsinomiss[1,7]
mydfB$VP3S[aaa]<-totalsinomiss[2,6]
mydfB$VP3M[aaa]<-totalsinomiss[2,7]
mydfB$VP4S[aaa]<-totalsinomiss[3,6]
mydfB$VP4M[aaa]<-totalsinomiss[3,7]
mydfB$VP1S[aaa]<-totalsinomiss[4,6]
mydfB$VP1M[aaa]<-totalsinomiss[4,7]
mydfB$VP2S[aaa]<-totalsinomiss[5,6]
mydfB$VP2M[aaa]<-totalsinomiss[5,7]
mydfB$ElicasiS[aaa]<-totalsinomiss[6,6]
mydfB$ElicasiM[aaa]<-totalsinomiss[6,7]
mydfB$VPgS[aaa]<-totalsinomiss[7,6]
mydfB$VPgM[aaa]<-totalsinomiss[7,7]
mydfB$treCRdRpS[aaa]<-totalsinomiss[8,6]
mydfB$treCRdRpM[aaa]<-totalsinomiss[8,7]
}

write.table(mydfB, file="C:/Documents/mutazioni_Sinmiss/allgenesDWVB.txt", quote=FALSE, row.names=FALSE, sep="\t")

#si calcola la media delle frequenze di mutazione sinonime e missenso per ciascun gene per il virus A, salvato in un file allgenesDWVB



#GRAFICI E TEST VIRUS B SINONIME PER TUTTI I GENI - 03/06

setwd("C:/Documents/png")
png("boxplot_sinonime_DWV-B.png",units="cm",res=600, width=15, height=15)
virus_B_sin_allgenes<-boxplot(mydfB$LpS,mydfB$VP3S,mydfB$VP4S,mydfB$VP1S,mydfB$VP2S,mydfB$ElicasiS,mydfB$VPgS,mydfB$treCRdRpS, main="Boxplot mutazioni sinonime DWV-B", names=c("LpS","VP3S","VP4S","VP1S","VP2S","ElicasiS","VPgS","3C-RdRpS"), ylab="Gene diversity",las=2,ylim=c(0,0.05))
virus_B_sin_allgenes
dev.off()


#GRAFICI E TEST VIRUS B MISSENSO PER TUTTI I GENI - 03/06

png("boxplot_missenso_DWV-B.png",units="cm",res=600, width=15, height=15)
virus_B_miss_allgenes<-boxplot(mydfB$LpM,mydfB$VP3M,mydfB$VP4M,mydfB$VP1M,mydfB$VP2M,mydfB$ElicasiM,mydfB$VPgM,mydfB$treCRdRpM, main="Boxplot mutazioni missenso DWV-B", names=c("LpM","VP3M","VP4M","VP1M","VP2M","ElicasiM","VPgM","3C-RdRpM"), ylab="Gene diversity",las=2, ylim=c(0,0.03))
virus_B_miss_allgenes
dev.off()



#TABELLA MEDIANA PER OGNI GENE PER DWV-A e DWV-B SINONIME E MISSENSO - 08/06

setwd("C:/Documents/mediana")

virus_A_sin_allgenes$stats[3,]
virus_A_miss_allgenes$stats[3,]
virus_B_sin_allgenes$stats[3,]
virus_B_miss_allgenes$stats[3,]
mediana_table<-data.frame(CAMPIONE= c("DWV-A MEDIANA SINONIME", "DWV-A MEDIANA MISSENSO", "DWV-B MEDIANA SINONIME", "DWV-B MEDIANA MISSENSO"), Lp=NA, VP3=NA,VP4=NA,VP1=NA,VP2=NA,Elicasi=NA,VPg=NA,treCRdRp=NA,stringsAsFactors = FALSE)

mediana_table[1,2:9]<-virus_A_sin_allgenes$stats[3,1:8]
mediana_table [2,2:9]<-virus_A_miss_allgenes$stats[3,1:8]
mediana_table [3,2:9]<-virus_B_sin_allgenes$stats[3,1:8]
mediana_table [4,2:9]<-virus_B_miss_allgenes$stats[3,1:8]
mediana_table
#calcola la mediana in una tabella per i 4 casi dei boxplot finali

write.table(mediana_table, file="mediana_table.txt", quote=FALSE,row.names=FALSE,sep="\t")





#GRAFICI E TEST STATISTICI VIRUS A VS VIRUS B, mutazionii sinonime - 19/06

setwd("C:/Documents/mutazioni_Sinmiss")

allgenesDWVA<-fread("allgenesDWVA.txt",data.table=FALSE)
allgenesDWVB<-fread("allgenesDWVB.txt",data.table=FALSE)

W_LpS_AvsB<-wilcox.test(allgenesDWVA$LpS,allgenesDWVB$LpS)
#W = 2146, p-value = 0.0001173
W_LpS_BvsA<-wilcox.test(allgenesDWVB$LpS,allgenesDWVA$LpS)
#W = 862, p-value = 0.0001173, quindi A>B
W_VP3S_AvsB<-wilcox.test(allgenesDWVA$VP3S,allgenesDWVB$VP3S)
#W = 1277, p-value = 0.1728
W_VP3S_BvsA<-wilcox.test(allgenesDWVB$VP3S,allgenesDWVA$VP3S)
#W = 1731, p-value = 0.1728, quindi B>A
W_VP4S_AvsB<-wilcox.test(allgenesDWVA$VP4S,allgenesDWVB$VP4S)
#W = 1927.5, p-value = 0.00764
W_VP4S_BvsA<-wilcox.test(allgenesDWVB$VP4S,allgenesDWVA$VP4S)
#W = 1080.5, p-value = 0.00764, quindi A>B
W_VP1S_AvsB<-wilcox.test(allgenesDWVA$VP1S,allgenesDWVB$VP1S)
#W = 1309.5, p-value = 0.2446
W_VP1S_BvsA<-wilcox.test(allgenesDWVB$VP1S,allgenesDWVA$VP1S)
#W = 1698.5, p-value = 0.2446, quindi B>A
W_VP2S_AvsB<-wilcox.test(allgenesDWVA$VP2S,allgenesDWVB$VP2S)
#W = 867.5, p-value = 0.0001262
W_VP2S_BvsA<-wilcox.test(allgenesDWVB$VP2S,allgenesDWVA$VP2S)
#W = 2140.5, p-value = 0.0001262, quindi B>A
W_ElicasiS_AvsB<-wilcox.test(allgenesDWVA$ElicasiS,allgenesDWVB$ElicasiS)
#W = 1043.5, p-value = 0.005798
W_ElicasiS_BvsA<-wilcox.test(allgenesDWVB$ElicasiS,allgenesDWVA$ElicasiS)
#W = 1964.5, p-value = 0.005798, quindi B>A
W_VPgS_AvsB<-wilcox.test(allgenesDWVA$VPgS,allgenesDWVB$VPgS)
#W = 1135, p-value = 0.02415
W_VPgS_BvsA<-wilcox.test(allgenesDWVB$VPgS,allgenesDWVA$VPgS)
#W = 1873, p-value = 0.02415, quindi B>A
W_treCRdRpS_AvsB<-wilcox.test(allgenesDWVA$treCRdRpS,allgenesDWVB$treCRdRpS)
#W = 997.5, p-value = 0.002469
W_treCRdRpS_BvsA<-wilcox.test(allgenesDWVB$treCRdRpS,allgenesDWVA$treCRdRpS)
#W = 2010.5, p-value = 0.002469, quindi B>A


#GRAFICI E TEST STATISTICI VIRUS A VS VIRUS B, mutazioni missenso - 19/06


W_LpM_AvsB<-wilcox.test(allgenesDWVA$LpM,allgenesDWVB$LpM)
#W = 2639, p-value = 8.955e-12
W_LpM_BvsA<-wilcox.test(allgenesDWVB$LpM,allgenesDWVA$LpM)
#W = 369, p-value = 8.955e-12, quindi A>B
W_VP3M_AvsB<-wilcox.test(allgenesDWVA$VP3M,allgenesDWVB$VP3M)
#W = 2621.5, p-value = 3.849e-13
W_VP3M_BvsA<-wilcox.test(allgenesDWVB$VP3M,allgenesDWVA$VP3M)
#W = 386.5, p-value = 3.849e-13, quindi A>B
W_VP4M_AvsB<-wilcox.test(allgenesDWVA$VP4M,allgenesDWVB$VP4M)
#W = 1527.5, p-value = 0.4016
W_VP4M_BvsA<-wilcox.test(allgenesDWVB$VP4M,allgenesDWVA$VP4M)
#W = 1480.5, p-value = 0.4016, quindi A>B
W_VP1M_AvsB<-wilcox.test(allgenesDWVA$VP1M,allgenesDWVB$VP1M)
#W = 2469, p-value = 3.361e-09
W_VP1M_BvsA<-wilcox.test(allgenesDWVB$VP1M,allgenesDWVA$VP1M)
#W = 539, p-value = 3.361e-09, quindi A>B
W_VP2M_AvsB<-wilcox.test(allgenesDWVA$VP2M,allgenesDWVB$VP2M)
#W = 2504, p-value = 5.19e-11
W_VP2M_BvsA<-wilcox.test(allgenesDWVB$VP2M,allgenesDWVA$VP2M)
#W = 504, p-value = 5.19e-11, quindi A>B
W_ElicasiM_AvsB<-wilcox.test(allgenesDWVA$ElicasiM,allgenesDWVB$ElicasiM)
#W = 2338, p-value = 2.627e-07
W_ElicasiM_BvsA<-wilcox.test(allgenesDWVB$ElicasiM,allgenesDWVA$ElicasiM)
#W = 670, p-value = 2.627e-07, quindi A>B
W_VPgM_AvsB<-wilcox.test(allgenesDWVA$VPgM,allgenesDWVB$VPgM)
#W = 2637, p-value = 4.494e-14
W_VPgM_BvsA<-wilcox.test(allgenesDWVB$VPgM,allgenesDWVA$VPgM)
#W = 371, p-value = 4.494e-14, quindi A>B
W_treCRdRpM_AvsB<-wilcox.test(allgenesDWVA$treCRdRpM,allgenesDWVB$treCRdRpM)
#W = 2371.5, p-value = 1.539e-07
W_treCRdRpM_BvsA<-wilcox.test(allgenesDWVB$treCRdRpM,allgenesDWVA$treCRdRpM)
#W = 636.5, p-value = 1.539e-07, quindi A>B


#GRAFICI E TEST STATISTICI VIRUS A, MUTAZIONI MISSENSO LP VS ALTRI GENI - 20/06

W_LpM_VS_VP3M<-wilcox.test(allgenesDWVA$LpM,allgenesDWVA$VP3M)
# p-value = 1.282e-13
W_LpM_VS_VP4M<-wilcox.test(allgenesDWVA$LpM,allgenesDWVA$VP4M)
#p-value < 2.2e-16
W_LpM_VS_VP1M<-wilcox.test(allgenesDWVA$LpM,allgenesDWVA$VP1M)
#p-value = 1.573e-13
W_LpM_VS_VP2M<-wilcox.test(allgenesDWVA$LpM,allgenesDWVA$VP2M)
#p-value = 1.099e-13
W_LpM_VS_ElicasiM<-wilcox.test(allgenesDWVA$LpM,allgenesDWVA$ElicasiM)
#p-value = 1.572e-13
W_LpM_VS_VPgM<-wilcox.test(allgenesDWVA$LpM,allgenesDWVA$VPgM)
#p-value = 5.002e-14
W_LpM_VS_treCRdRpM<-wilcox.test(allgenesDWVA$LpM,allgenesDWVA$treCRdRpM)
#p-value = 8.258e-13

#GRAFICI E TEST STATISTICI VIRUS B, MUTAZIONI MISSENSO LP VS ALTRI GENI - 20/06

WB_LpM_VS_VP3M<-wilcox.test(allgenesDWVB$LpM,allgenesDWVB$VP3M)
# p-value = 2.259e-11
WB_LpM_VS_VP4M<-wilcox.test(allgenesDWVB$LpM,allgenesDWVB$VP4M)
#p-value = 5.808e-13
WB_LpM_VS_VP1M<-wilcox.test(allgenesDWVB$LpM,allgenesDWVB$VP1M)
#p-value = 1.101e-06
WB_LpM_VS_VP2M<-wilcox.test(allgenesDWVB$LpM,allgenesDWVB$VP2M)
#p-value = 7.711e-11
WB_LpM_VS_ElicasiM<-wilcox.test(allgenesDWVB$LpM,allgenesDWVB$ElicasiM)
#p-value = 3.545e-06
WB_LpM_VS_VPgM<-wilcox.test(allgenesDWVB$LpM,allgenesDWVB$VPgM)
#p-value = 1.639e-07
WB_LpM_VS_treCRdRpM<-wilcox.test(allgenesDWVB$LpM,allgenesDWVB$treCRdRpM)
#p-value = 0.002453