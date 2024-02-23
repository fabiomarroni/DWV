#fare un loop dei file csv aggiornati, e unirli per posizione (merge(fileone, filetwo, by ="POS", all=TRUE, sort=TRUE))

library(dplyr)
library(data.table)

# Specifica il percorso della cartella contenente i file 
cartella <- "/Users/michelegarbellotti/Desktop/R_genetic/data_DWV_gene_finestre"

# Ottieni l'elenco dei file csv nella cartella
files <- list.files(path = cartella, pattern = "DWV-A", full.names = TRUE)

# Inizializza il data frame risultante con il primo file
result_df <- read.csv(files[1])

# Itera sui file rimanenti e uniscili al data frame risultante
for (file in files[-1]) {
  # Leggi il file CSV
  csv_data <- read.csv(file)
  
  # Unisci i dataframes mantenendo l'ordine delle colonne
  result_df <- merge(result_df, csv_data, by = "POS", all=TRUE, sort=TRUE)
}

# Colonnes da rimuovere
colonne_da_rimuovere <- c("FILTER", "ID", "FREQ_REF", "GENE_DIV", "FREQ_REF_MEDIA")

# Rimuovi le colonne indesiderate
result_df <- result_df %>%
  select(-one_of(colonne_da_rimuovere))

# Salva il data frame risultante in un unico file
output_file <- paste0(cartella, "_unito.csv")
write.csv(result_df, file = output_file, row.names = FALSE)

# Seleziona le colonne che iniziano con "FREQ_ALT_MEDIA"
colonne_freq_alt_media <- grep("^FREQ_ALT_MEDIA", names(result_df), value = TRUE)

# Calcola la media per ogni colonna selezionata
media_totale_freq_alt_media <- colMeans(result_df[, colonne_freq_alt_media], na.rm = TRUE)

# Calcola la media del vettore di valori medi
freq_alt <- mean(media_totale_freq_alt_media)

# Ottieni l'elenco dei file csv nella cartella
files1 <- list.files(path = cartella, pattern = "DWV-B", full.names = TRUE)

# Inizializza il data frame risultante con il primo file
result_dfB <- read.csv(files1[1])

# Itera sui file rimanenti e uniscili al data frame risultante
for (file in files1[-1]) {
  # Leggi il file CSV
  csv_dataB <- read.csv(file)
  
  # Unisci i dataframes mantenendo l'ordine delle colonne
  result_dfB <- merge(result_dfB, csv_dataB, by = "POS", all=TRUE, sort=TRUE)
}

# Colonnes da rimuovere
colonne_da_rimuovereB <- c("FILTER", "ID", "FREQ_REF", "GENE_DIV", "FREQ_REF_MEDIA")

# Rimuovi le colonne indesiderate
result_dfB <- result_dfB %>%
  select(-one_of(colonne_da_rimuovereB))

# Salva il data frame risultante in un unico file
output_file1 <- paste0(cartella, "_unitoB.csv")
write.csv(result_dfB, file = output_file1, row.names = FALSE)

# Seleziona le colonne che iniziano con "FREQ_ALT_MEDIA"
colonne_freq_alt_mediaB <- grep("^FREQ_ALT_MEDIA", names(result_df), value = TRUE)

# Calcola la media per ogni colonna selezionata
media_totale_freq_alt_mediaB <- colMeans(result_dfB[, colonne_freq_alt_mediaB], na.rm = TRUE)

# Calcola la media del vettore di valori medi
freq_altB <- mean(media_totale_freq_alt_mediaB)
