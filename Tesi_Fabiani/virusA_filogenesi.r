#FILOGENESI VIRUS A

library(ape)
library(RColorBrewer)
setwd("C:/Documents/fasta/")
filoA<-read.dna(file="totalA.fasta",format = "fasta",as.matrix = FALSE,as.character = FALSE)
#non mi esegue il forceDnatolower per avere le basi in maiuscolo
filoA

setwd("C:/Documents/etichette_virus_nuove")
getwd()
namesA<-fread(file="etichette api tesi_MF.txt")
namesA
namesA$"NOMI SAMPLE DWV-A"<-gsub(".VCF","",namesA$"NOMI SAMPLE DWV-A")
namesA
names(filoA)[!is.na(match(names(filoA),namesA$"NOMI SAMPLE DWV-A"))]= na.omit(namesA$"NUOVE ETICHETTE DWV-A"[match(names(filoA),namesA$"NOMI SAMPLE DWV-A")])
names(filoA)
#modfica i nomi con le nuove etichette

setwd("C:/Documents/fasta/")
class(filoA)
#la classe è Dnabin
write.dna(filoA, "filoA.fasta", format = "fasta", nbcol = -1, colsep = "")
#salva il fasta
distA<-dist.dna(filoA)
#calcola la distanza tra le sequenze di DNA
length(distA)
#calcola la lunghezza
alberoA<-bionj(distA)
#si tratta di un algoritmo per costruire l'albero 
alberoA
write.tree(alberoA,"../Newick/alberoA.nwk")
#salva l'albero filogenetico
plot(alberoA,cex=.5,"u")
#mi da come output l'albero filogenetico

#VERIFICA NJ
library(stats)
setwd("C:/Documents/png/")
x <- as.vector(distA)
y <- as.vector(as.dist(cophenetic(alberoA)))
png("njA.png")
plot(x, y, xlab="Distanze coppie originali", ylab="Distanze coppie sull'albero", main="VERIFICA CORRETTEZZA NJ", pch=20, cex=3, xlim=c(0,0.025),ylim=c(0,0.025))
abline(lm(y~x), col="red")
dev.off()
#salvo il grafico dell'algoritmo in un file png
cor(x,y)^2
#verifica se il dataset è compatibile con questo algoritmo

#HEATMAP
mycolA<-brewer.pal(9, "YlOrRd")
#comando che permette di scegliere i colori da usare poi nella heatmapù
setwd("C:/Documents/png")
png("heatmapA.png",units="cm",res=600, width=15, height=12)
heatmap(as.matrix(distA),col=rev(mycolA),margins=c(7,9),cexRow=0.25, cexCol=0.25, symm=T)
#col=rev serve a invertire i colori in modo che scuro significhi vicino
dev.off()
#salva la heat map all'interno di un pdf

#ALGORITMO UPGMA
library(phangorn)
alberoAupgma<- upgma(distA,method="average")
plot(alberoAupgma,cex=.5,"u")

#VERIFICA UPGMA
library(stats)
png("UPGMA_A.png")
treeA <- as.phylo(upgma(distA,method="average"))
x <- as.vector(distA)
y <- as.vector(as.dist(cophenetic(alberoAupgma)))
plot(x, y, xlab="Distanze coppie originali", ylab="Distanze coppie sull'albero", main="VERIFICA CORRETTEZZA UPGMA", pch=20,cex=3,xlim=c(0,0.025),ylim=c(0,0.025))
abline(lm(y~x), col="red")
dev.off()
cor(x,y)^2
#verifica se l'algoritmo UPGMA è corretto in questo caso

#DENDROGRAMMA A CLUSTER
alberoAupgma3<- hclust(distA,method="average")
plot(alberoBupgma3,cex=.5)
#non spiega come verificarne la correlazione
