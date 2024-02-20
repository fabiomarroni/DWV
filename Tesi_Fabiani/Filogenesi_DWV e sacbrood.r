#ANALISI FILOGENETICA REFERENCE A E B CON SACBROOD 1 E 2


setwd("C:/Documents/fasta/")
#la directory dev'essere impostata dentro alla cartella fasta
library(seqinr)

myseq<-read.alignment(file="DWVA_B_SAC1e2.aln",format="clustal",forceToLower=F)
#il file clustal deriva dal sito clustalW e dall'unione delle ref di virus A,B e i due SAC
myseqmat<-as.matrix.alignment(myseq)

myseqmat<-myseqmat[,848:9932]
#evidenzia le posizioni da 848 a 9932 ovvero dove devo tagliare la sequenza

myseqstring<-apply(myseqmat,1,paste,sep="",collapse="")
write.fasta(as.list(myseqstring),row.names(myseqmat),file.out="DWV_sacbrood_alignment.fasta",as.string=T)
#salvo il fasta contenete le 4 sequenze tagliate+

#remove.gap(file="DWV_sacbrood_alignment",file.out="DWV_sacbrood_alig")
#per il momento non è necessario

library(ape)
library(RColorBrewer)
filoDWVsac<-read.dna(file="DWV_sacbrood_alignment.fasta",format = "fasta",as.matrix = FALSE,as.character = FALSE)
#non mi esegue il forceDnatolower per avere le basi in maiuscolo
filoDWVsac
class(filoDWVsac)
#la classe è Dnabin
write.dna(filoDWVsac, "filoDWVsac.fasta", format = "fasta", nbcol = -1, colsep = "")
#salva il fasta
distDWVsac<-dist.dna(filoDWVsac)
#calcola la distanza tra le sequenze di DNA
length(distDWVsac)
#calcola la lunghezza
alberoDWVsac<-bionj(distDWVsac)
#si tratta di un algoritmo per costruire l'albero 
alberoDWVsac
write.tree(alberoDWVsac,"../Newick/alberoDWVsac.nwk")
#salva l'albero filogenetico
plot(alberoDWVsac,cex=.5,"u")
#mi da come output l'albero filogenetico

#VERIFICA NJ
library(stats)
x <- as.vector(distDWVsac)
y <- as.vector(as.dist(cophenetic(alberoDWVsac)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is NJ appropriate?", pch=20, cex=3)
abline(lm(y~x), col="red")
cor(x,y)^2
#verifica se il dataset è compatibile con questo algoritmo

#HEATMAP
mycolDWVsac<-brewer.pal(9, "YlOrRd")
#comando che permette di scegliere i colori da usare poi nella heatmapù
setwd("C:/Documents/pdf")
pdf("heatmapDWVsac.pdf", paper="a4")
heatmap(as.matrix(distDWVsac),col=rev(mycolDWVsac))
#col=rev serve a invertire i colori in modo che scuro significhi vicino
dev.off()
#salva la heat map all'interno di un pdf


#ANALISI FILOGENETICA DWV-A CON SACBROOD 1 E 2 E REFERENCE DI B

library(seqinr)
setwd("C:/Documents/fasta/")
#la directory dev'essere impostata dentro alla cartella fasta

myseq<-read.alignment(file="tutti_i_DWVA_refB_SAC1e2.aln",format="clustal",forceToLower=F)
#il file clustal deriva dal sito clustalW e dall'unione di tutti i virus A, la ref di B e i due SAC
myseqmat<-as.matrix.alignment(myseq)

myseqmat<-myseqmat[,848:9932]
#evidenzia le posizioni da 848 a 9932 ovvero dove devo tagliare la sequenza

myseqstring<-apply(myseqmat,1,paste,sep="",collapse="")
write.fasta(as.list(myseqstring),row.names(myseqmat),file.out="DWVA_refB_sacbrood_alignment.fasta",as.string=T)
#salvo il fasta contenete le sequenze tagliate

library(ape)
library(RColorBrewer)
filoDWVAsac<-read.dna(file="DWVA_refB_sacbrood_alignment.fasta",format = "fasta",as.matrix = FALSE,as.character = FALSE)
#non mi esegue il forceDnatolower per avere le basi in maiuscolo
filoDWVAsac
class(filoDWVAsac)
#la classe è Dnabin
write.dna(filoDWVAsac, "filoDWVAsac.fasta", format = "fasta", nbcol = -1, colsep = "")
#salva il fasta
distDWVAsac<-dist.dna(filoDWVAsac)
#calcola la distanza tra le sequenze di DNA
length(distDWVAsac)
#calcola la lunghezza
alberoDWVAsac<-bionj(distDWVAsac)
#si tratta di un algoritmo per costruire l'albero 
alberoDWVAsac
write.tree(alberoDWVAsac,"../Newick/alberoDWVAsac.nwk")
#salva l'albero filogenetico
plot(alberoDWVAsac,cex=.5,"u")

#VERIFICA NJ
library(stats)
x <- as.vector(distDWVAsac)
y <- as.vector(as.dist(cophenetic(alberoDWVAsac)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is NJ appropriate?", pch=20, cex=3)
abline(lm(y~x), col="red")
cor(x,y)^2
#verifica se il dataset è compatibile con questo algoritmo

#HEATMAP
mycolDWVAsac<-brewer.pal(9, "YlOrRd")
#comando che permette di scegliere i colori da usare poi nella heatmapù
setwd("C:/Documents/pdf")
pdf("heatmapDWVAsac.pdf", paper="a4")
heatmap(as.matrix(distDWVAsac),col=rev(mycolDWVAsac))
#col=rev serve a invertire i colori in modo che scuro significhi vicino
dev.off()
#salva la heat map all'interno di un pdf




#ANALISI FILOGENETICA DWV-B CON SACBROOD 1 E 2 E REFERENCE DI A

library(seqinr)
setwd("C:/Documents/fasta/")
#la directory dev'essere impostata dentro alla cartella fasta

myseq<-read.alignment(file="tutti_i_DWVB_refA_SAC1e2.aln",format="clustal",forceToLower=F)
#il file clustal deriva dal sito clustalW e dall'unione di tutti i virus A, la ref di B e i due SAC
myseqmat<-as.matrix.alignment(myseq)

myseqmat<-myseqmat[,848:9932]
#evidenzia le posizioni da 848 a 9932 ovvero dove devo tagliare la sequenza

myseqstring<-apply(myseqmat,1,paste,sep="",collapse="")
write.fasta(as.list(myseqstring),row.names(myseqmat),file.out="DWVB_refA_sacbrood_alignment.fasta",as.string=T)
#salvo il fasta contenete le sequenze tagliate

library(ape)
library(RColorBrewer)
filoDWVBsac<-read.dna(file="DWVB_refA_sacbrood_alignment.fasta",format = "fasta",as.matrix = FALSE,as.character = FALSE)
#non mi esegue il forceDnatolower per avere le basi in maiuscolo
filoDWVBsac
class(filoDWVBsac)
#la classe è Dnabin
write.dna(filoDWVBsac, "filoDWVBsac.fasta", format = "fasta", nbcol = -1, colsep = "")
#salva il fasta
distDWVBsac<-dist.dna(filoDWVBsac)
#calcola la distanza tra le sequenze di DNA
length(distDWVBsac)
#calcola la lunghezza
alberoDWVBsac<-bionj(distDWVBsac)
#si tratta di un algoritmo per costruire l'albero 
alberoDWVBsac
write.tree(alberoDWVBsac,"../Newick/alberoDWVBsac.nwk")
#salva l'albero filogenetico
plot(alberoDWVBsac,cex=.5,"u")

#VERIFICA NJ
library(stats)
x <- as.vector(distDWVBsac)
y <- as.vector(as.dist(cophenetic(alberoDWVBsac)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is NJ appropriate?", pch=20, cex=3)
abline(lm(y~x), col="red")
cor(x,y)^2
#verifica se il dataset è compatibile con questo algoritmo

#HEATMAP
mycolDWVBsac<-brewer.pal(9, "YlOrRd")
#comando che permette di scegliere i colori da usare poi nella heatmapù
setwd("C:/Documents/pdf")
pdf("heatmapDWVBsac.pdf", paper="a4")
heatmap(as.matrix(distDWVBsac),col=rev(mycolDWVBsac),symm=T)
#col=rev serve a invertire i colori in modo che scuro significhi vicino
dev.off()
#salva la heat map all'interno di un pdf


#ANALISI FILOGENETICA DWVA, DWVB CON SACBROOD 1 E 2 E REFERENCE DI A e B - 23/05

library(seqinr)
setwd("C:/Documents/fasta/")
#la directory dev'essere impostata dentro alla cartella fasta

myseq<-read.alignment(file="tutti_DWVAeB_refAeB_SAC1e2.aln",format="clustal",forceToLower=F)
#il file clustal deriva dal sito clustalW e dall'unione di tutti i virus A, la ref di B e i due SAC
myseqmat<-as.matrix.alignment(myseq)

myseqmat<-myseqmat[,848:9932]
#evidenzia le posizioni da 848 a 9932 ovvero dove devo tagliare la sequenza



#ANALISI FILOGENETICA REF A E B REF SAC, 3 CAMPIONI PER DWV-A E 3 CAMPIONI PER DWV-B

library(seqinr)
setwd("C:/Documents/fasta/")
#la directory dev'essere impostata dentro alla cartella fasta

myseq<-read.alignment(file="refAeB_refSAC1e2_6campioniAeB.aln",format="clustal",forceToLower=F)
#il file clustal deriva dal sito clustalW e dall'unione di tutti i virus A, la ref di B e i due SAC
myseq$nam[1]<-"DWV-A_reference"
myseq$nam[2]<-"DWV-A_ISA_2011_3"
myseq$nam[3]<-"DWV-A_ISNb_2009_2"
myseq$nam[4]<-"DWV-A_USNb_2009_1"
myseq$nam[5]<-"DWV-B_reference"
myseq$nam[6]<-"DWV-B_ISA_2011_1"
myseq$nam[7]<-"DWV-B_ISNb_2009_8"
myseq$nam[8]<-"DWV-B_USNb_2009_1"
myseq$nam[9]<-"SAC-1_reference"
myseq$nam[10]<-"SAC-2_reference"
#modifico i nomi in base alle etichette nuove

myseqmat<-as.matrix.alignment(myseq)

myseqmat<-myseqmat[,848:9932]
#evidenzia le posizioni da 848 a 9932 ovvero dove devo tagliare la sequenza

myseqstring<-apply(myseqmat,1,paste,sep="",collapse="")
write.fasta(as.list(myseqstring),row.names(myseqmat),file.out="refAeB_refSAC1e2_6campioniAeB.fasta",as.string=T)
#salvo il fasta contenete le sequenze tagliate

library(ape)
library(RColorBrewer)
filoDWVABsac<-read.dna(file="refAeB_refSAC1e2_6campioniAeB.fasta",format = "fasta",as.matrix = FALSE,as.character = FALSE)
#non mi esegue il forceDnatolower per avere le basi in maiuscolo
filoDWVABsac
class(filoDWVABsac)
#la classe è Dnabin
write.dna(filoDWVABsac, "refAeB_refSAC1e2_6campioniAeB_alignment.fasta", format = "fasta", nbcol = -1, colsep = "")
#salva il fasta
distDWVABsac<-dist.dna(filoDWVABsac)
#calcola la distanza tra le sequenze di DNA
length(distDWVABsac)
#calcola la lunghezza
alberoDWVABsac<-bionj(distDWVABsac)
#si tratta di un algoritmo per costruire l'albero 
alberoDWVABsac
write.tree(alberoDWVABsac,"../Newick/alberoDWVAeB_6campioni_sac.nwk")
#salva l'albero filogenetico
plot(alberoDWVABsac,cex=.5)

#VERIFICA NJ
library(stats)
x <- as.vector(distDWVABsac)
y <- as.vector(as.dist(cophenetic(alberoDWVABsac)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is NJ appropriate?", pch=20, cex=3)
abline(lm(y~x), col="red")
cor(x,y)^2
#verifica se il dataset è compatibile con questo algoritmo

#HEATMAP
mycolDWVABsac<-brewer.pal(9, "YlOrRd")
#comando che permette di scegliere i colori da usare poi nella heatmapù
setwd("C:/Documents/png")
png("heatmapDWVAeB_6campioni_sac.png",units="cm",res=600, width=15, height=12)
heatmap(as.matrix(distDWVABsac),col=rev(mycolDWVABsac),margins=c(7,9),cexRow=0.7, cexCol=0.7,symm=T)
#col=rev serve a invertire i colori in modo che scuro significhi vicino
dev.off()
#salva la heat map all'interno di un pdf