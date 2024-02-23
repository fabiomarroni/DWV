rm(list=ls())
library(vcfR)
library(data.table)
library(dplyr)
library(VariantAnnotation)

# Definisci la funzione che esegue il tuo procedimento
procedimento <- function(myDWV) {
  
  DWV_1 <- data.frame(vcf_data@fix)
  DWV_40 <- DWV_1[grep("40.000", DWV_1$QUAL, fixed = TRUE), ]
  
  myDWV <- (DWV_40$INFO)
  
  #divido le strighe in base agli MM creando liste di valori MM es 123,44,555  12,34
  split1 <- strsplit(myDWV, ";")
  parts <- unlist(split1)
  mm_part <- grep("^MM=", parts, value = TRUE) #risultato tutti i MMs=xx
  mm_value <- sub("MM=", "", mm_part) #risultato solo i valori xx
  mm_value <- data.frame(mm_value)
  #mm_max <- sapply(strsplit(mm_value, ","), function(x) max(as.numeric(x))) prende il valore di mm più grande per ogni lista di valori
  #mm_max <- data.frame(mm_max)
  
  #divido le stringhe in base alla virgola creando liste di basi es a,t,g  c,g
  secondDWV <- DWV_40$ALT
  ALT_part <- data.frame(secondDWV)
  
  #unisco i due gruppi di liste creando un dataframe a due colonne 
  DWV_Data <- data.frame(Colonna1 = mm_value, Colonna2 = ALT_part)
  colnames(DWV_Data) <- c("MM_MAX", "ALT_MAX")
  
  # Funzione per ottenere il massimo valore da una riga
  trova_massimo <- function(row) {
    
    basi_alt <- unlist(strsplit(as.character(row["ALT_MAX"]), ","))
    valori_mm <- as.numeric(unlist(strsplit(as.character(row["MM_MAX"]), ",")))
    
    indice_massimo <- which.max(valori_mm)
    valore_massimo <- valori_mm[indice_massimo]
    base_corrispondente <- basi_alt[indice_massimo]
    
    return(data.frame(ValoreMassimoMM = valore_massimo, BaseCorrispondente = base_corrispondente))
  }
  
  # Applica la funzione a ciascuna riga del dataframe
  risultati <- do.call(rbind, apply(DWV_Data, 1, trova_massimo))
  
  # Aggiungi i risultati al dataframe originale
  DWV_Data$MM_MAX <- risultati[, 1]
  DWV_Data$ALT_MAX <- risultati[, 2]
  
  #Creo valori MM
  mm_max <- grep("", DWV_Data$MM_MAX, value = TRUE)
  mm_numeric <- as.numeric(mm_max)
  mm_numeric[is.na(mm_numeric)] <- 0 
  is.numeric(mm_numeric)
  
  # Cerco i valori che iniziano con DP che fanno parte delle stringhe che contengono anche MMsum
  dp_part1 <- grep("MMsum", myDWV, value = TRUE) #trovo le stringhe che contengono MMsum e di conseguenza hanno anche DP (dato che è presente in ogni riga) 
  split1 <- unlist(strsplit(dp_part1, ";")) #le divido in base al punto e virgola e le rendo non liste
  dp_part2 <- grep("DP=", split1, value = TRUE) #ricavo solo le stringhe con i valori di DP 
  dp_numeric1 <- as.numeric(sub("DP=", "", dp_part2)) #determino solo i numeri di DP 
  dp_numeric1[is.na(dp_numeric1)] <- 0 
  
  # calcolo le frequenze degli alleli
  freq_ALT <- round( signif( (mm_numeric/dp_numeric1), digits = 2), digits = 2) 
  freq_REF <- 1-freq_ALT 
  
  # Creare un altro dataframe con le nuove colonne delle frequenze
  DWV_freq <- data.frame(Colonna1 = freq_ALT, Colonna2 = freq_REF)
  colnames(DWV_freq) <- c("FREQ_ALT", "FREQ_REF")
  
  # Considera le righe che soddisfano una condizione (la colonna INFO dve contenere MM)
  #DWV_MM <- DWV_40[grep("MM", DWV_40$INFO, fixed = TRUE), ]
  
  # Calcolare la gene diversity
  gene_div <- round( signif( (2*freq_ALT*freq_REF), digits = 2), digits = 2)
  gene_div[gene_div ==  0] <- NA
  DWV_gene <- data.frame(gene_div)
  colnames(DWV_gene) <- c("GENE_DIV")
  
  # Unisco i dataframe delle frequenze e quelli della gene diversity
  DWV_freq <- cbind(DWV_freq, DWV_gene)
  
  calcola_medie <- function(DWV_freq, variabile, lunghezza_finestra = 100, passo = 10) {
    lunghezza_dataframe <- nrow(DWV_freq)
    
    # Inizializza un vettore per le medie
    medie <- numeric(length = lunghezza_dataframe)
    
    n_finestra <- numeric(length = lunghezza_dataframe)
    
    # Controlla se il dataframe ha meno di 100 elementi
    if (lunghezza_dataframe < lunghezza_finestra) {
      cat("Il dataframe ha meno di 100 elementi. Restituendo la media delle frequenze disponibili.")
      return(list(valori_medie = mean(DWV_freq[[variabile]], na.rm = TRUE), valori_n = 1))
      
    }
    
    # Calcola la media delle frequenze in finestre scorrevoli
    for (i in seq(1, lunghezza_dataframe - lunghezza_finestra + 1, by = passo)) {
      finestra <- DWV_freq[[variabile]][i:(i + lunghezza_finestra - 1)]
      media_finestra <- mean(finestra, na.rm = TRUE)
      medie[i:(i + lunghezza_finestra - 1)] <- rep(media_finestra, length.out = lunghezza_finestra)
      n_finestra[i:(i + lunghezza_finestra - 1)] <- rep(DWV_40$POS[i], length.out = lunghezza_finestra)
      
    }
    
    return(list(valori_medie = medie, valori_n = n_finestra))
  }
  
  # Utilizzo della funzione per calcolare le medie
  FREQ_ALT_MEDIA <- calcola_medie(DWV_freq, "FREQ_ALT")$valori_medie
  FREQ_REF_MEDIA <- calcola_medie(DWV_freq, "FREQ_REF")$valori_medie
  GENE_DIV_MEDIA <- calcola_medie(DWV_freq, "GENE_DIV")$valori_medie
  
  N_FINESTRA <- calcola_medie(DWV_freq, "FREQ_ALT")$valori_n
  
  # Creazione dei dataframe
  FREQ_ALT_MEDIA <- signif(data.frame(FREQ_ALT_MEDIA), digits = 2)
  FREQ_REF_MEDIA <- signif(data.frame(FREQ_REF_MEDIA), digits = 2)
  GENE_DIV_MEDIA <- signif(data.frame(GENE_DIV_MEDIA), digits = 2)
  
  # Assegnazione dei nomi alle colonne
  colnames(FREQ_ALT_MEDIA) <- c("FREQ_ALT_MEDIA")
  colnames(FREQ_REF_MEDIA) <- c("FREQ_REF_MEDIA")
  colnames(GENE_DIV_MEDIA) <- c("GENE_DIV_MEDIA")
  
# Unire i dataframe
DWV_new <- cbind(DWV_40, DWV_freq, FREQ_ALT_MEDIA, FREQ_REF_MEDIA, GENE_DIV_MEDIA, N_FINESTRA)

# sostituisco le basi alternative che hanno MM più grande
DWV_new[, c("ALT")] <- DWV_Data$ALT_MAX
  
return(DWV_new)
}

# Specifica il percorso della cartella contenente i file VCF
cartella <- "/Users/michelegarbellotti/Desktop/R_genetic/03_vir_SNP"

# Ottieni l'elenco dei file VCF nella cartella
files <- list.files(path = cartella, pattern = "DWV", full.names = TRUE)

# Verifico se la lista creata continene i file DWV
print(files)

# Loop attraverso ogni file
for (file in files) {
  
  # Leggi il file VCF
  vcf_data <- read.vcfR(file)
  
  # Applica la funzione al campo desiderato del dataframe vcf_data
  result <- procedimento(vcf_data)
  
  # salvali in un nuovo file
  SNP_aggiornato <- paste0(tools::file_path_sans_ext(file), "_aggiornato.csv")
  
  write.csv(result, SNP_aggiornato, row.names = FALSE) 
  
  result_rid <- result[!duplicated(result$N_FINESTRA), ]
  
  write.csv(result_rid, gsub(".csv", "_rid.csv", SNP_aggiornato), row.names = FALSE)
}


