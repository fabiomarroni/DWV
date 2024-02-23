# Contenuto  script R principale
source("/Users/michelegarbellotti/Desktop/DWV-FINALE.R")
library(ggplot2)

tajima.DWV_A_valori <- as.numeric(unlist(tajima.DWV_A))
tajima.DWV_B_valori <- as.numeric(unlist(tajima.DWV_B))

finposmean_valori <-as.numeric(unlist(finposmean))
finposmeanB_valori <-as.numeric(unlist(finposmeanB))

titvmean_valori <- as.numeric(unlist(titvmean))
titvmeanB_valori <- as.numeric(unlist(titvmeanB))

tab_valoriA <- data.frame(POS = finposmean_valori, TAJIMA = tajima.DWV_A_valori)
tab_valoriB <- data.frame(POS = finposmeanB_valori, TAJIMA = tajima.DWV_B_valori)

tab2_valoriA <- data.frame(POS = finposmean_valori, RAPP = titvmean_valori)
tab2_valoriB <- data.frame(POS = finposmeanB_valori, RAPP = titvmeanB_valori)

mean_rapp_A <- summary(tab2_valoriA$RAPP)
mean_rapp_B <- summary(tab2_valoriB$RAPP)
mean_Tajima_A <- mean(tab_valoriA$TAJIMA)
mean_Tajima_B <- mean(tab_valoriB$TAJIMA)

#CONFRONTO TRA CEPPI A E B

# Aggiungi una colonna per identificare il ceppo
tab_valoriA$Ceppo <- "A"
tab_valoriB$Ceppo <- "B"

png("boxplot_tajima.png", width = 800, height = 600)
boxplot(tab_valoriA$TAJIMA, tab_valoriB$TAJIMA,
        xlab = "Ceppo ",
        ylab = "Tajima",
        main = "Boxplot di confronto tra i valori di Tajima'D dei ceppi A e B",
        col= c("blue", "red"),
        names = c("A", "B"))
points(x = jitter(rep(1, length(regioni_significative$TAJIMA)), factor = 0.1), y = regioni_significative$TAJIMA, col = "black", pch = 1)
points(x = jitter(rep(2, length(regioni_significative$TAJIMA)), factor = 0.1), y =  regioni_significative$TAJIMA, col = "black", pch = 1)
# Chiudi il dispositivo grafico
dev.off()


#IDENTIFICAZIONE DI REGIONI SIGNIFICATIVE

# Trova le posizioni con i valori di Tajima più alti e più bassi per il ceppo A
max_tajima_A <- tab_valoriA[order(tab_valoriA$TAJIMA, decreasing = TRUE), ][1:4, ]
min_tajima_A <- tab_valoriA[order(tab_valoriA$TAJIMA), ][1:4, ]

# Trova le posizioni con i valori di Tajima più alti e più bassi per il ceppo B
max_tajima_B <- tab_valoriB[order(tab_valoriB$TAJIMA, decreasing = TRUE), ][1:4, ]
min_tajima_B <- tab_valoriB[order(tab_valoriB$TAJIMA), ][1:4, ]

# Crea un dataframe con i risultati
regioni_significative <- rbind(
  data.frame(Ceppo = "A", Tipo = "Massimo", max_tajima_A),
  data.frame(Ceppo = "A", Tipo = "Minimo", min_tajima_A),
  data.frame(Ceppo = "B", Tipo = "Massimo", max_tajima_B),
  data.frame(Ceppo = "B", Tipo = "Minimo", min_tajima_B)
)

#Grafico Tajima ceppi A e B

#processo per salvarlo come foto
png("grafico_tajima.png", width = 800, height = 600)

# Crea il grafico per DWV_A
plot(finposmean, as.numeric(unlist(tajima.DWV_A)), col = "blue", type = "l", lty = 1, lwd = 2, xlab = "posizione nucleotidica", ylab = "Valore di Tajima's D", ylim = c(-1.5, 1.5))
# Aggiungi il grafico per DWV_B
lines(finposmeanB, as.numeric(unlist(tajima.DWV_B)), col = "red", type = "l", lty = 1, lwd = 2)
# Aggiungi la legenda
legend("topright", legend = c("DWV_A", "DWV_B"), col = c("blue", "red"), lty = c(1, 1), lwd = 2)
# Imposta i valori per l'asse x (ogni 500)
ticks_x <- seq(0, max_x, by = 500)
# Aggiungi valori sull'asse x
axis(side = 1, at = ticks_x, labels = ticks_x)
# Aggiungi un titolo
title(main = "Confronto dei valori di Tajima tra DWV_A e DWV_B")
points(regioni_significative$POS, regioni_significative$TAJIMA, col = "black", pch = 16, cex = 0.5)

# Chiudi il dispositivo grafico
dev.off()


#Grafico Ti/Tv ceppi A e B

# Apri il dispositivo grafico PNG
png("grafico_Ti_Tv.png", width = 800, height = 600)
# Crea il grafico per il ceppo A
plot(finposmean, as.numeric(unlist(titvmean[,3])), col = "blue", type = "l", lty = 1, lwd = 2, xlab = "posizione nucleotidica", ylab = "Ti/TV", ylim = c(0, 15))
# Aggiungi il grafico per il ceppo B
lines(finposmeanB, as.numeric(unlist(titvmeanB[,3])), col = "red", type = "l", lty = 1, lwd = 2)
# Aggiungi la legenda
# Imposta i valori per l'asse x (ogni 500)
ticks_x <- seq(0, max_x, by = 500)
# Aggiungi valori sull'asse x
axis(side = 1, at = ticks_x, labels = ticks_x)
legend("topright", legend = c("Ti/Tv ceppo A", "Ti/Tv ceppo B"), col = c("black", "red"), lty = c(1, 1), lwd = 2)
# Aggiungi un titolo
title(main = "Confronto dei rapporti Ti/Tv tra DWV_A e DWV_B")
# Chiudi il dispositivo grafico
dev.off()

