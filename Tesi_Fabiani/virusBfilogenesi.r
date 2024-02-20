#ANALISI FILOGENETICA VIRUS B, ALGORITMI, HEATMAP

#ALGORITMO BIONJ

library(ape)
library(RColorBrewer)
library (data.table)
setwd("C:/Documents/fasta/")
filoB<-read.dna(file="totalB.fasta",format = "fasta",as.matrix = FALSE,as.character = FALSE)
#non mi esegue il forceDnatolower per avere le basi in maiuscolo
filoB

setwd("C:/Documents/etichette_virus_nuove")
namesB<-fread(file="etichette api tesi_MF.txt")
namesB
namesB$"NOMI SAMPLE DWV-B"<-gsub(".VCF","",namesB$"NOMI SAMPLE DWV-B")
namesB
names(filoB)[!is.na(match(names(filoB),namesB$"NOMI SAMPLE DWV-B"))]= na.omit(namesB$"NUOVE ETICHETTE DWV-B"[match(names(filoB),namesB$"NOMI SAMPLE DWV-B")])
names(filoB)
#modifica le etichette con i nuovi nomi

setwd("C:/Documents/fasta/")
class(filoB)
#la classe è Dnabin
write.dna(filoB, "filoB.fasta", format = "fasta", nbcol = -1, colsep = "")
#salva il fasta
distB<-dist.dna(filoB)
#calcola la distanza tra le sequenze di DNA
length(distB)
#calcola la lunghezza
alberoB<-bionj(distB)
#si tratta di un algoritmo per costruire l'albero 
alberoB
write.tree(alberoB,"../Newick/alberoB.nwk")
#salva l'albero filogenetico
plot(alberoB,cex=.5,"u")
#mi da come output l'albero filogenetico con "u" che serve a disegnare un albero unrooted

#VERIFICA NJ
library(stats)
setwd("C:/Documents/png/")
x <- as.vector(distB)
y <- as.vector(as.dist(cophenetic(alberoB)))
png("njB.png")
plot(x, y, xlab="Distanze coppie originali", ylab="Distanze coppie sull'albero", main="VERIFICA CORRETTEZZA NJ", pch=20, cex=3, xlim=c(0,0.04),ylim=c(0,0.04))
abline(lm(y~x), col="red")
dev.off()
#salvo il grafico dell'algoritmo in un file png
cor(x,y)^2
#verifica se il dataset è compatibile con questo algoritmo

#HEATMAP
mycolB<-brewer.pal(9, "YlOrRd")
#comando che permette di scegliere i colori da usare poi nella heatmap
setwd("C:/Documents/png")
png("heatmapB.png",units="cm",res=600, width=15, height=12)
heatmap(as.matrix(distB),col=rev(mycolB),margins=c(7,9),cexRow=0.25, cexCol=0.25, symm=T)
#col=rev serve a invertire i colori in modo che scuro significhi vicino
dev.off()
#salva la heat map in un png

#ALGORITMO UPGMA
library(phangorn)
alberoBupgma<- upgma(distB,method="average")
plot(alberoBupgma,cex=.5,"u")

#VERIFICA UPGMA
library(stats)
png("UPGMA_B.png")
treeB <- as.phylo(upgma(distB,method="average"))
x <- as.vector(distB)
y <- as.vector(as.dist(cophenetic(alberoBupgma)))
plot(x, y, xlab="Distanze coppie originali", ylab="Distanze coppie sull'albero", main="VERIFICA CORRETTEZZA UPGMA", pch=20,cex=3, xlim=c(0,0.04),ylim=c(0,0.04))
abline(lm(y~x), col="red")
dev.off()
cor(x,y)^2
#verifica se l'algoritmo UPGMA è corretto in questo caso

#DENDROGRAMMA A CLUSTER
alberoBupgma3<- hclust(distB,method="average")
plot(alberoBupgma3,cex=.5)
#non spiega come verificarne la correlazione

#ALTRE HEATMAP

install.packages("gplots")
library("gplots")
df<- scale(distB)
heatmap.2(df, scale = "none", col = rev(bluered(100)),trace = "none", density.info = "none")
#si tratta di una heatmap migliorata

