rm(list=ls())
library(dplyr)
library(data.table)
library(VariantAnnotation)
library(tidyverse)

# Specifica il percorso della cartella contenente i file VCF
cartella <- "/Users/michelegarbellotti/Desktop/R_genetic/data_DWV_freq_qual_gene"

files_csv <- list.files(path = cartella, pattern = "\\.csv$", full.names = TRUE)

# Inizializza un vettore per contare il numero di righe per ogni file CSV
num_righe <- numeric(length(files_csv))
gene_div_sum <- numeric(length(files_csv))

# Loop attraverso i file CSV e controlla il numero di righe
for (i in seq_along(files_csv)) {
  data <- read.csv(files_csv[i])
  num_righe[i] <- nrow(data)
  gene_div_sum[i] <- mean(data$GENE_DIV, na.rm = TRUE)
}

# Crea un dataframe con i risultati
n_SNP <- data.frame(File = basename(files_csv), NumRows = num_righe, GeneDivSum = gene_div_sum)
n_SNP[is.na(n_SNP)] <- 0

# raccogli i file che iniziano con un determinato codice
n_SNP_DWV_A_78 <- n_SNP %>%
  filter(grepl("DWV-A_78", File)) 
n_SNP_DWV_B_78 <- n_SNP %>%
  filter(grepl("DWV-B_78", File))
n_SNP_DWV_B_78 <- n_SNP_DWV_B_78 %>%
  slice(-4, -5, -7, -13, -19)

n_SNP_DWV_B_75 <- n_SNP %>%
  filter(grepl("DWV-B_75", File)) 

n_SNP_DWV_B_57 <- n_SNP %>%
  filter(grepl("DWV-B_57", File)) 

n_SNP_DWV_A_Sample <- n_SNP %>%
  filter(grepl("DWV-A_Sample", File)) 

n_SNP_DWV_A_ape <- n_SNP %>% #11
  filter(grepl("DWV-A_ape", File)) 

n_SNP_DWV_A_APE <- n_SNP %>%
  filter(grepl("DWV-A_APE", File)) 


file_path <- "/Users/michelegarbellotti/Desktop/R_genetic/Virus_A_B.txt"

#leggere il file in formato txt
Virus_A_B <- read.table(file_path, header = TRUE)

# raccogli i file che iniziano con un determinato codice
Virus_A_B_57 <- Virus_A_B %>%
  filter(grepl("^57", File) & DWVgp1 > 0)
Virus_A_B_57 <- Virus_A_B_57 %>% 
  slice(-n())
Virus_A_B_75 <- Virus_A_B %>%
  filter(grepl("^75", File) & DWVgp1 > 1) 

Virus_A_B_78 <- Virus_A_B %>%
  filter(grepl("^78", File)) 
Virus_A_B_78 <- Virus_A_B_78 %>% 
  slice(-4, -5, -7, -13, -19, -34)

Virus_A_B_ape <- Virus_A_B %>%
  filter(grepl("^ape", File))  

Virus_A_B_APE <- Virus_A_B %>%
  filter(grepl("^APE", File))  

Virus_A_B_Sample <- Virus_A_B %>%
  filter(grepl("^Sample", File)) 
Virus_A_B_Sample <- Virus_A_B_Sample %>% 
  slice(-3)

#unione in due dataframe in base al ceppo del virus
n_SNP_DWV_A <- bind_rows(n_SNP_DWV_A_78, n_SNP_DWV_A_ape, n_SNP_DWV_A_APE, n_SNP_DWV_A_Sample)
n_SNP_DWV_B <- bind_rows(n_SNP_DWV_B_57, n_SNP_DWV_B_75, n_SNP_DWV_B_78)

#unione delle reads in due dataframe in base al ceppo del virus
Virus_A <- bind_rows(Virus_A_B_78, Virus_A_B_ape, Virus_A_B_APE, Virus_A_B_Sample)
Virus_B <- bind_rows(Virus_A_B_57, Virus_A_B_75, Virus_A_B_78)

#coefficente correlazione tra le due misure
coef_cor_DWV_B_57 <- cor(n_SNP_DWV_B_57$NumRows, Virus_A_B_57$VDV1_gp1)
coef_cor_DWV_B_75 <- cor(n_SNP_DWV_B_75$NumRows, Virus_A_B_75$VDV1_gp1)
coef_cor_DWV_B_78 <- cor(n_SNP_DWV_B_78$NumRows, Virus_A_B_78$VDV1_gp1)
coef_cor_DWV_A_78 <- cor(n_SNP_DWV_A_78$NumRows, Virus_A_B_78$DWVgp1)
coef_cor_DWV_A_Sample <- cor(n_SNP_DWV_A_Sample$NumRows, Virus_A_B_Sample$DWVgp1)
coef_cor_DWV_A_ape <- cor(n_SNP_DWV_A_ape$NumRows, Virus_A_B_ape$DWVgp1)
coef_cor_DWV_A_APE <- cor(n_SNP_DWV_A_APE$NumRows, Virus_A_B_APE$DWVgp1)

coef_cor_DWV_B2_57 <- cor(n_SNP_DWV_B_57$GeneDivSum, Virus_A_B_57$VDV1_gp1)
coef_cor_DWV_B2_75 <- cor(n_SNP_DWV_B_75$GeneDivSum, Virus_A_B_75$VDV1_gp1)
coef_cor_DWV_B2_78 <- cor(n_SNP_DWV_B_78$GeneDivSum, Virus_A_B_78$VDV1_gp1)
coef_cor_DWV_A2_78 <- cor(n_SNP_DWV_A_78$GeneDivSum, Virus_A_B_78$DWVgp1)
coef_cor_DWV_A2_Sample <- cor(n_SNP_DWV_A_Sample$GeneDivSum, Virus_A_B_Sample$DWVgp1)
coef_cor_DWV_A2_ape <- cor(n_SNP_DWV_A_ape$GeneDivSum, Virus_A_B_ape$DWVgp1)
coef_cor_DWV_A2_APE <- cor(n_SNP_DWV_A_APE$GeneDivSum, Virus_A_B_APE$DWVgp1)

#calcolo coefficienti di correlazione 
coef_cor_DWV_A <- cor(n_SNP_DWV_A$NumRows, Virus_A$DWVgp1)
coef_cor_DWV_B <- cor(n_SNP_DWV_B$NumRows, Virus_B$VDV1_gp1)

coef_cor_DWV_A2 <- cor(n_SNP_DWV_A$GeneDivSum, Virus_A$DWVgp1)
coef_cor_DWV_B2 <- cor(n_SNP_DWV_B$GeneDivSum, Virus_B$VDV1_gp1)

#dataframe con numero SNP e reads
df_DWV_A <- data.frame(File = n_SNP_DWV_A$File, SNP = n_SNP_DWV_A$NumRows, Reads_Virus = Virus_A$DWVgp1)
df_DWV_B <- data.frame(File = n_SNP_DWV_B$File, SNP = n_SNP_DWV_B$NumRows, Reads_Virus = Virus_B$VDV1_gp1)

#dataframe con gene diversity e reads
df_DWV_A2 <- data.frame(File = n_SNP_DWV_A$File, GENE = n_SNP_DWV_A$GeneDivSum, Reads_Virus = Virus_A$DWVgp1)
df_DWV_B2 <- data.frame(File = n_SNP_DWV_B$File, GENE = n_SNP_DWV_B$GeneDivSum, Reads_Virus = Virus_B$VDV1_gp1)

library(ggplot2)
library(ggrepel)

theme_custom <- theme_minimal() +
  theme(
    # Modifica delle dimensioni del testo per migliorare la leggibilità
    text = element_text(size = 12),
    # Impostazione dei colori
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgray"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"  # Rimozione della legenda
  )

#Grafico DWV-A Diversità genetica
graficoA2 <- ggplot(df_DWV_A2, aes(x = GENE, y = Reads_Virus, label = File)) +
  geom_point(color = "blue") +  # Punti blu
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linea di tendenza blu
  geom_text(aes(x = max(GENE), y = max(Reads_Virus), 
                label = paste("Correlazione:", round(coef_cor_DWV_A2, 2))),
            hjust = 1.2, vjust = 1.2, color = "black") +  # Testo sulla retta di tendenza blu
  labs(title = "Relazione tra Diversità nucleotidica e Reads di Virus",
       x = "Valore Diversità nucleotidica",
       y = "Numero di Reads di Virus") +
  scale_y_continuous(labels = scales::comma_format()) +  # Formattazione delle etichette sull'asse y
  theme_custom
ggsave("grafico_DWV_A_gene.png", plot = graficoA2, width = 8, height = 6, units = "in")

#Grafico DWV_A snp
graficoA <- ggplot(df_DWV_A, aes(x = SNP, y = Reads_Virus, label = File)) +
  geom_point(color = "blue") +  # Punti blu
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linea di tendenza blu
  geom_text(aes(x = max(SNP), y = max(Reads_Virus), 
                label = paste("Correlazione:", round(coef_cor_DWV_A, 2))),
            hjust = 1.2, vjust = 1.2, color = "black") +  # Testo sulla retta di tendenza blu
  labs(title = "Relazione tra SNP e Reads di Virus",
       x = "Numero di SNP",
       y = "Numero di Reads di Virus") +
  scale_y_continuous(labels = scales::comma_format()) +  # Formattazione delle etichette sull'asse y
  theme_custom

ggsave("grafico_DWV_A_snp.png", plot = graficoA, width = 8, height = 6, units = "in")

#Grafico DWV_B snp
graficoB <- ggplot(df_DWV_B, aes(x = SNP, y = Reads_Virus, label = File)) +
  geom_point(color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Aggiungi la linea di tendenza
  geom_text(aes(x = max(SNP), y = max(Reads_Virus), 
            label = paste("Correlazione:", round(coef_cor_DWV_B, 2))), 
            hjust = 1.2, vjust = 1.2, color = "black") +  # Aggiungi il testo sulla retta di tendenza
  labs(title = "Relazione tra SNP e Reads di Virus",
       x = "Numero di SNP",
       y = "Numero di Reads di Virus") +
  scale_y_continuous(labels = scales::comma_format()) +
  theme_custom

ggsave("grafico_DWV_B.png", plot = graficoB, width = 8, height = 6, units = "in")

#Grafico DWV-B diversità nucleotidica
graficoB2 <- ggplot(df_DWV_B2, aes(x = GENE, y = Reads_Virus, label = File)) +
  geom_point(color = "red") +  # Punti blu
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linea di tendenza blu
  geom_text(aes(x = max(GENE), y = max(Reads_Virus), 
                label = paste("Correlazione:", round(coef_cor_DWV_B2, 2))),
            hjust = 1.2, vjust = 1.2, color = "black") +  # Testo sulla retta di tendenza blu
  labs(title = "Relazione tra Diversità nucleotidica e Reads di Virus",
       x = "Valore Diversità nucleotidica",
       y = "Numero di Reads di Virus") +
  scale_y_continuous(labels = scales::comma_format()) +  # Formattazione delle etichette sull'asse y
  theme_custom

ggsave("grafico_DWV_B_gene.png", plot = graficoB2, width = 8, height = 6, units = "in")





