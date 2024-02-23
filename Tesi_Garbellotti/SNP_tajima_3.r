##fare loop per fineste di 100 per calcolo medie frequenze alt, fare grafico tajima stile (type plot l)
rm(list=ls())
install.packages("dplyr")
library(dplyr)

#file_path <- result_df
#file_pathB <- result_dfB

file_path <- "/Users/michelegarbellotti/Desktop/R_genetic/data_DWV_gene_finestre_unito.csv"
file_pathB <- "/Users/michelegarbellotti/Desktop/R_genetic/data_DWV_gene_finestre_unitoB.csv"

datafame_merge <- read.csv(file_path)
datafame_mergeB <- read.csv(file_pathB)

#Dato che facciamo la media nel loop, adesso prendiamo la GD di ogni SNP, e non la media
datafame_merge <- select(datafame_merge, !contains("MEDIA"))
datafame_mergeB <- select(datafame_mergeB, !contains("MEDIA"))

#Calcolo approssimativo della GD a livello di popolazione (media delle GD individuali) 
datafame_merge$APPROX_GD<-rowMeans(select(datafame_merge,starts_with("GENE_DIV")),na.rm=T)
datafame_mergeB$APPROX_GD<-rowMeans(select(datafame_mergeB,starts_with("GENE_DIV")),na.rm=T)

#Le posizioni in cui APPROX_GD è NA sono prive di SNP, per cui le buttiamo (forse non è importante, ma meglio farlo)
datafame_merge<-datafame_merge[!is.na(datafame_merge$APPROX_GD),]
datafame_mergeB<-datafame_mergeB[!is.na(datafame_mergeB$APPROX_GD),]

#Stessa cosa per quando si ha GD=0
datafame_merge<-datafame_merge[datafame_merge$APPROX_GD>0,]
datafame_mergeB<-datafame_mergeB[datafame_mergeB$APPROX_GD>0,]

# Calcola il numero totale di righe nel DataFrame
total_rows <- nrow(datafame_merge)
total_rowsB <- nrow(datafame_mergeB)

# Calcolo delle basi reference e alternative
tref<-data.frame(select(datafame_merge,starts_with("REF")))
trefB<-data.frame(select(datafame_mergeB,starts_with("REF")))

tref[is.na(tref)]<-""
trefB[is.na(trefB)]<-""

talt<-data.frame(select(datafame_merge,starts_with("ALT")))
taltB<-data.frame(select(datafame_mergeB,starts_with("ALT")))

talt[is.na(talt)]<-""
taltB[is.na(taltB)]<-""

datafame_merge$ALLREF<-substr(apply(apply(tref,1,unlist),2,paste0,collapse=""),1,1)
datafame_mergeB$ALLREF<-substr(apply(apply(trefB,1,unlist),2,paste0,collapse=""),1,1)

#Qua si seleziona come allele alternativo il primo che si incontra.
#Il prossimo passo potrebbe essere vedere nella stringa quante volte si trova ciascuna base e poi usare la base più frequente
datafame_merge$ALLALT<-substr(apply(apply(talt,1,unlist),2,paste0,collapse=""),1,1)
datafame_mergeB$ALLALT<-substr(apply(apply(taltB,1,unlist),2,paste0,collapse=""),1,1)

#Calcolo delle proprorzioni per trovare Ti e TV e loro rapporto
calcola_proporzioni <- function(REF, ALT) {
  # Creazione di un vettore contenente coppie di basi REF e ALT
  coppie_basi <- paste0(REF, ALT)
  
  # Definizione di coppie di transizione e trasversione
  transizioni <- c("AG", "GA", "CT", "TC")
  trasversioni <- c("AC", "CA", "AT", "TA", "GC", "CG", "GT", "TG")
  
  # Calcolo del numero di transizioni e trasversioni
  num_transizioni <- sum(coppie_basi %in% transizioni)
  num_trasversioni <- sum(coppie_basi %in% trasversioni)
  
  # Calcolo delle proporzioni
  proporzione_transizioni <- num_transizioni / length(coppie_basi)
  proporzione_trasversioni <- num_trasversioni / length(coppie_basi)
  
  # Creazione di un data.frame con i risultati
  risultati <- data.frame(
    Proporzione_Transizioni_Ti = proporzione_transizioni,
    Proporzione_Trasversioni_Tv = proporzione_trasversioni,
    rapp_TiTv <- proporzione_transizioni/proporzione_trasversioni)
  
  numeri <- data.frame(
    num_transizioni <- sum(coppie_basi %in% transizioni),
    num_trasversioni <- sum(coppie_basi %in% trasversioni))
  
  return(list(risultati = risultati, numeri =numeri))
}

# Dimensione della finestra e dimensione del passo
window_size <- 100
step_size <- 15

# Inizializza un vettore per memorizzare le medie
newmeans<-finposmean <- titvmean <- nnmean<- numeric()
newmeansB<-finposmeanB <- titvmeanB <- nnmeanB<- numeric()

# Loop attraverso il DataFrame
for (i in seq(1, total_rows - window_size + 1, by = step_size)) {
  
  # Seleziona le prime 100 righe nella finestra corrente
  current_window <- datafame_merge[i:(i + window_size - 1), ]
 
  # Seleziona le colonne che iniziano con "gene_div"
  gene_div_columns <- select(current_window, starts_with("GENE_DIV"))
  
  # Calcola la media delle colonne selezionate
  current_mean <- colMeans(gene_div_columns, na.rm = TRUE)
  approx_gd_mean<-mean(current_window$APPROX_GD,na.rm=T)
  pos_mean<-mean(current_window$POS,na.rm=T)
  tt_mean<-calcola_proporzioni(current_window$ALLREF,current_window$ALLALT)$risultati
  nn_mean<-calcola_proporzioni(current_window$ALLREF,current_window$ALLALT)$numeri
  
  
  # Aggiungi la media al vettore
  newmeans <- rbind(newmeans,approx_gd_mean)
  finposmean<-rbind(finposmean,pos_mean)
  titvmean<-rbind(titvmean,tt_mean)
  nnmean<-rbind(nnmean, nn_mean)
}

for (i in seq(1, total_rowsB - window_size + 1, by = step_size)) {
  
  # Seleziona le prime 100 righe nella finestra corrente
  current_windowB <- datafame_mergeB[i:(i + window_size - 1), ]
  
  # Seleziona le colonne che iniziano con "freq_alt_media"
  freq_alt_columnsB <- select(current_windowB, starts_with("GENE_DIV"))
  
  # Calcola la media delle colonne selezionate
  current_meanB <- colMeans(freq_alt_columnsB, na.rm = TRUE)
  approx_gd_meanB<-mean(current_windowB$APPROX_GD,na.rm=T)
  pos_meanB<-mean(current_windowB$POS,na.rm=T)
  tt_meanB<-calcola_proporzioni(current_windowB$ALLREF,current_windowB$ALLALT)$risultati
  nn_meanB<-calcola_proporzioni(current_windowB$ALLREF,current_windowB$ALLALT)$numeri
  
  # Aggiungi la media al vettore
  newmeansB <- rbind(newmeansB,approx_gd_meanB)
  finposmeanB<-rbind(finposmeanB,pos_meanB)
  titvmeanB<-rbind(titvmeanB,tt_meanB)
  nnmeanB<-rbind(nnmeanB, nn_meanB)
}

### CALCOLO DEL TAJIMA ####
cartella <- "/Users/michelegarbellotti/Desktop/R_genetic/data_DWV_gene_finestre"
##S, che è quanti SNP ci sono nella finestra, 
S <- window_size

##nn quanti cromosomi (per noi praticamente quanti file di input del virus DWVA o DWVB, 
files1 <- list.files(path = cartella, pattern = "DWV-A", full.names = TRUE)
files2 <- list.files(path = cartella, pattern = "DWV-B", full.names = TRUE)
nn <- length(files1)
nnB <- length(files2)

##e poi tajima.gd in cui metteremo la gene diversity, che tu sai già calcolare, ma dobbiamo calcolarla non per ogni singolo individuo, bensì per tutta la popolazione.
##Potremmo cominciare usando come gene diversity la media della gene diversity di una finestra (che puoi calcolare adattando la funzione della frequenza delle medie), e poi (dato che avrai un valore per ogni individuo), dare come tajima.gd la media tra tutti gli individui.

# Funzione per il calcolo di Tajima's D
tajima.d<-function(S, nn, tajima.gd)
{
  a1<-sum(1/seq(1:nn))
  a2<-sum(1/(seq(1:nn)^2))
  tw<-S/a1
  b1<-(nn+1)/(3*(nn-1))
  b2<-2*(nn^2+nn+3)/(9*nn*(nn-1))
  c1<-b1-(1/a1)
  e1<-c1/a1
  c2<-b2-((nn+2)/(a1*nn))+a2/a1^2
  e2<-c2/(a1^2+a2)
  tajima.d<-(tajima.gd-tw)/sqrt(e1*S+e2*S*(S-1))
  tajima.d
}

tajima.DWV_A <- data.frame(tajima.d(S, nn, S*newmeans))
tajima.DWV_B <- data.frame(tajima.d(S, nnB, S*newmeansB))




