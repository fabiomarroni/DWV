#23/02 - VIRUS B SOSTITUZIONE SNP E SALVATAGGIO FASTA

setwd("C:/Documents/SNP")
#imposta la directory nella cartella SNP
library("Biostrings")
library("seqinr")
library(Rsubread)
library(data.table)
myfilesB<-list.files(path=".", pattern=".VCF", all.files=FALSE, full.names=FALSE)
#crea una lista con tutti i file nella cartella della directory, in questo caso VCF


myfilesB<-myfilesB[grep("DWV-B",myfilesB)]
#devo selezionare solo le tabelle del virus B
highcovB<-fread("../ABcount/ABcountvirusB.txt")$Sample
highcovB<-gsub("ReadsPerGene.out.tab","",highcovB)
newfilesB<-rep(NA,length(highcovB))
for(aaa in 1:length(highcovB))
{
newfilesB[aaa]<-myfilesB[grep(highcovB[aaa],myfilesB)]
}
#il loop serve a tenere solo i file con alto coverage >5000 dal file di testo creato con ABCcount
myfilesB<-newfilesB
write(myfilesB,"myfilesB.txt")
#salva la lista come file di testo
for(aaa in 1:length(myfilesB))
{
a<-myfilesB[aaa]
#a considera scorre tra tutti i files e costituisce un nome univoco per tutti
print(a)
VCF<-fread(a,skip=9)
#carica tutti i files come tabella che chiama VCF
if (nrow(VCF)<1) next
#serve a saltare i file VCF vuoti
VCF$INFO<-gsub("MM=","",VCF$INFO)
#elimina la scritta MM= dalla colona INFO
VCF$cov1<-NA
VCF$cov2<-NA
VCF$cov3<-NA
#servono a creare le colonne cov1...con valore NA
VCF$INFO<-gsub("INDEL","",VCF$INFO)
#serve a eliminare la scritta INDEL (inutile)
VCF$allcov<-unlist(lapply(strsplit(VCF$INFO,";"),"[",3))
#inserisco gli elementi che vengono 3 volte dopo il; della colonna INFO nella colonna allcov
prova<-"1;2;3"
strsplit(prova,";")
#il comando strsplit separa gli elementi
unlist(strsplit(prova,";"))
#forma una lista separando gli elementi
lapply(strsplit(prova,";"),"[",2)
#fornisce come output il secondo numero dopo;
lapply(strsplit(prova,";"),"[",3)
VCF$cov1<-unlist(lapply(strsplit(VCF$allcov,","),"[",1))
#mette nella colonna cov1 il primo valore della colonna allcov
VCF$cov1<-gsub("SR=","",VCF$cov1)
#elimino la scritta SR= dalla colonna cov1
VCF$allcov<-gsub("SR=","",VCF$allcov)
VCF$cov2<-unlist(lapply(strsplit(VCF$allcov,","),"[",2))
#inserisco il secondo valore della colonna allcov nella colonna cov2
VCF$cov3<-unlist(lapply(strsplit(VCF$allcov,","),"[",3))
#inserisco il terzo valore della colonna allcov nella colonna cov3
VCF$allcov<-NULL
#mi elimina la colonna indicata
VCF$DP<-unlist(lapply(strsplit(VCF$INFO,";"),"[",1))
#prende il primo valore della colonna INFO e lo riporta nella colonna DP
VCF$DP<-gsub("DP=","",VCF$DP)
#elimina la scritta DP= dalla colonna DP
VCF$Freq.<-NA
#crea la colonna frequenza di mutazione
VCF$SNP1<-NA
VCF$SNP2<-NA
VCF$SNP3<-NA
#ho creato le tre colonne SNP1,2,3
VCF$SNP1<-unlist(lapply(strsplit(VCF$ALT,","),"[",1))
VCF$SNP2<-unlist(lapply(strsplit(VCF$ALT,","),"[",2))
VCF$SNP3<-unlist(lapply(strsplit(VCF$ALT,","),"[",3))
#ho assegnato alle colonne SNP i rispettivi valori della colonna ALT
VCF$cov1=as.numeric(VCF$cov1)
#trasforma i caratteri testuali "" in numerici
sum(VCF$cov1)
#fa la somma dei valori della colonna cov1
VCF$DP=as.numeric(VCF$DP)
#trasforma i caratteri testuali in numerici della colonna DP
sum(VCF$DP)
#fa la somma dei valori della colonna cov1
VCF$DP<-unlist(lapply(strsplit(VCF$INFO,";"),"[",1))
# riinserisco i valori di DP dalla colonna INFO
VCF$DP<-gsub("DP=","",VCF$DP)
# elimino la scritta "DP=" dalla colonna DP
head(VCF)
VCF$QUAL>1
#seleziono tutto ciò che ha qualità maggiore di 1 in modo da tenere solo le SNP
VCF[VCF$QUAL>1,]
# tiene solo le righe contenenti le SNP ed elimina il resto
VCF=VCF[VCF$QUAL>1,]
nrow(VCF)
#vedo il numero di righe della tabella
VCF$DP=as.numeric(VCF$DP)
#trasformo la colonna DP in campo numerico 
head(VCF)
VCF$Freq.<-VCF$cov1/VCF$DP
#assegno alla colonna Freq. il risultato della divisione tra cov1 e DP
VCF
#mostra tutta la tabella
is.na(VCF$cov3)
VCF$cov3[is.na(VCF$cov3)]=0
#trova i valori NA della colonna cov3 e li trasforma in zero
VCF
VCF$cov2[is.na(VCF$cov2)]=0
#trova i valori NA della colonna cov2 e li trasforma in zero
VCF
VCF$Freq.<-VCF$cov1/VCF$DP*100
# ricalcolo le frequenze in %
VCF
VCF$cov2=as.numeric(VCF$cov2)
#trasforma i campi testuali della colonna cov2 in numerici
VCF$cov3=as.numeric(VCF$cov3)
#trasforma i campi testuali della colonna cov3 in numerici
VCF$Freq2<-VCF$cov2/VCF$DP*100
#mi calcola Freq2 in %
VCF$Freq3<-VCF$cov3/VCF$DP*100
#mi calcola Freq3 in %
VCF
VCF$newRef=VCF$REF
#aggiunge una colonna newRef che coincide con REF che sarà da modificare successivamente
VCF
for(bbb in 1:nrow(VCF))
{
if(VCF$Freq.[bbb]>=1 & VCF$Freq.[bbb]>VCF$Freq2[bbb] & VCF$Freq.[bbb]>VCF$Freq3[bbb]) VCF$newRef[bbb]=VCF$SNP1[bbb]
if(VCF$Freq2[bbb]>=1 & VCF$Freq2[bbb]>VCF$Freq.[bbb] & VCF$Freq2[bbb]>VCF$Freq3[bbb]) VCF$newRef[bbb]=VCF$SNP2[bbb]
if(VCF$Freq3[bbb]>=1 & VCF$Freq3[bbb]>VCF$Freq.[bbb] & VCF$Freq3[bbb]>VCF$Freq2[bbb]) VCF$newRef[bbb]=VCF$SNP3[bbb]
}
VCF
# quando il valore di Freq1 è maggiore di 1 e maggiore delle altre frequenze prende la base di SNP e la sostituisce nella colonna newREF, poi fa la stessa cosa con Freq2 e Freq3
VCF$newRef != VCF$REF
# idividua le righe dove la base azotata differisce tra la newRef e REF e mi da come risultato TRUE e FALSE
VCF[VCF$newRef != VCF$REF]
# mi da la tabella con solo le basi azotate che differiscono tra newRef e REF
Frequno=VCF[VCF$newRef != VCF$REF]
#creo una nuova tabella chiamata Frequno che rappresenta il dataframe per >1%
Frequno
VCF$newRef = VCF$REF
#reimposto la colona newRef uguale a REF per fare il calcolo con Freq>50%
VCF
for(ddd in 1:nrow(VCF))
{
if(VCF$Freq.[ddd]>50) VCF$newRef[ddd]=VCF$SNP1[ddd]
if(VCF$Freq2[ddd]>50) VCF$newRef[ddd]=VCF$SNP2[ddd]
if(VCF$Freq3[ddd]>50) VCF$newRef[ddd]=VCF$SNP3[ddd]
}
VCF
#questo è un loop: "aaa" è la variabile contatore numerica ed i valori che assume sono pari alle righe del dataframe, "in" cioè varia tra, 
#1:nrow(VCF) indica il numero di righe totali della tabella
#"if" se la frequenza è maggiore di 50%, "spazio" allora la colonna newref è uguale alla colonna SNP
VCF<-VCF[VCF$DP>=10,]
#Se non ci sono SNP con F>50% eDP>10 si salta al prossimo file
if(nrow(VCF)<1) next
#seleziona e tiene solo i valori dove DP è maggiore di 10 poichè dev'essere sufficientemente alta per avere una certa attendibilità
VCF$newRef != VCF$REF
VCF[VCF$newRef != VCF$REF]
Freqcinquanta=VCF[VCF$newRef != VCF$REF]
#creo una nuova tabella chiamata Freqcinquanta che rappresenta il dataframe per >50%
Freqcinquanta
write.table(Freqcinquanta,paste0("../Freqcinquanta/",a),quote=F,row.names=F)

#imposto la directory su fasta
#A<-read.fasta("virusA.fasta",forceDNAtolower = FALSE)

B<-read.fasta("virusB.fasta",forceDNAtolower = FALSE)
# è un comando del pacchetto seqinr che permette di leggere i file fasta e li leggo una volta per il virus A e una per il B

#visualizzo il file Fasta del virus A o B

#eseguo il loop per importare con il nome dei VCF automaticamentefor(aaa in 1:length(myfiles))

a<-myfilesB[aaa]
for(ccc in 1:nrow(Freqcinquanta))
{
if(B[[Freqcinquanta$"#CHROM"[ccc]]][Freqcinquanta$POS[ccc]]!=Freqcinquanta$REF[ccc] & (B[[Freqcinquanta$"#CHROM"[ccc]]][Freqcinquanta$POS[ccc]]%in% c("A","C","G","T"))) stop("Your reference is wrong") 
B[[Freqcinquanta$"#CHROM"[ccc]]][Freqcinquanta$POS[ccc]]=Freqcinquanta$newRef[ccc]
}

#permette di sostituire al file A o al file B la base delle SNP con Freq>50%
outfile<-gsub(".VCF",".fasta",a)
#write.fasta(A,names=names(A),file.out=paste0("../fasta/",outfile))
write.fasta(B,names=names(B),file.out=paste0("../fasta/",outfile))
#salva il fasta con lo stesso nome del VCF con le basi corrette
}

setwd("../fasta/")
DWVB<-list.files(pattern = "DWV-B")
#seleziona una lista con i fasta creati per il virus B

#DWVB<-myfiles[grep("DWV-B",myfiles)]

for(eee in 1:length(DWVB))
{
infasta<-read.fasta(file = DWVB[eee],forceDNAtolower=F)
seqname<-gsub(".fasta","",DWVB[eee])
if(eee==1)
{
totalB<-infasta
names(totalB)<-seqname
}
if(eee>1) 
{
names(infasta)<-seqname
totalB<-c(totalB,infasta)
}
}
write.fasta(sequences = totalB,names=names(totalB),file.out = "totalB.fasta")
#crea un unico file fasta chiamato totalB.fasta a partire da più file fasta
