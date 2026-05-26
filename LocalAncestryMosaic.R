#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

# Usage:
#   Rscript LocalAncestryMosaic.R <individual_id> <LC|HC> [base_directory]
#
# Examples:
#   Rscript LocalAncestryMosaic.R HG01976 LC /media/daniel/Espacio/NuevosBED
#   Rscript LocalAncestryMosaic.R HG01976 HC /media/daniel/Espacio/NuevosBED

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2 || length(args) > 3) {
  stop("Usage: Rscript LocalAncestryMosaic.R <individual_id> <LC|HC> [base_directory]")
}

individuo <- args[1]
set_datos <- toupper(args[2])
base_dir <- ifelse(length(args) == 3, args[3], "/media/daniel/Espacio/NuevosBED")

if(!set_datos %in% c("LC", "HC")) {
  stop("The second argument must be LC or HC")
}

etiqueta_set <- ifelse(set_datos == "LC", "lcWGS", "hcWGS")
carpeta_set <- paste0(set_datos, "0.9_GAP")

archivo_entrada <- file.path(
  base_dir,
  carpeta_set,
  paste0("InputG", individuo, ".tsv")
)

if(!file.exists(archivo_entrada)) {
  stop(paste("File does not exist:", archivo_entrada))
}

cat("Processing individual:", individuo, "\n")
cat("Dataset:", set_datos, "\n")
cat("Input file:", archivo_entrada, "\n")

datos <- read.csv(archivo_entrada, sep = "\t")

cromosomas <- paste0("chr", 1:22)

chr_sizes <- c(
  chr1=248956422,
  chr2=242193529,
  chr3=198295559,
  chr4=190214555,
  chr5=181538259,
  chr6=170805979,
  chr7=159345973,
  chr8=145138636,
  chr9=138394717,
  chr10=133797422,
  chr11=135086622,
  chr12=133275309,
  chr13=114364328,
  chr14=107043718,
  chr15=101991189,
  chr16=90338345,
  chr17=83257441,
  chr18=80373285,
  chr19=58617616,
  chr20=64444167,
  chr21=46709983,
  chr22=50818468
)

agregar_gaps_telomeros <- function(df, chr_size, chr) {
  if(nrow(df) == 0) {
    return(data.frame(
      Chr = chr,
      Inicio = 0,
      Fin = chr_size,
      colorhap1 = "gray",
      colorhap2 = "gray"
    ))
  }

  df <- df %>% arrange(Inicio)
  resultado <- data.frame()

  if(df$Inicio[1] > 0) {
    resultado <- bind_rows(
      resultado,
      data.frame(
        Chr = df$Chr[1],
        Inicio = 0,
        Fin = df$Inicio[1],
        colorhap1 = "gray",
        colorhap2 = "gray"
      )
    )
  }

  if(nrow(df) > 1) {
    for(i in 1:(nrow(df) - 1)) {
      resultado <- bind_rows(resultado, df[i, ])

      if(df$Fin[i] < df$Inicio[i + 1]) {
        resultado <- bind_rows(
          resultado,
          data.frame(
            Chr = df$Chr[i],
            Inicio = df$Fin[i],
            Fin = df$Inicio[i + 1],
            colorhap1 = "gray",
            colorhap2 = "gray"
          )
        )
      }
    }
  }

  resultado <- bind_rows(resultado, df[nrow(df), ])

  if(df$Fin[nrow(df)] < chr_size) {
    resultado <- bind_rows(
      resultado,
      data.frame(
        Chr = df$Chr[1],
        Inicio = df$Fin[nrow(df)],
        Fin = chr_size,
        colorhap1 = "gray",
        colorhap2 = "gray"
      )
    )
  }

  resultado
}

datos_completo <- data.frame()

for(chr in cromosomas) {
  tmp <- datos %>% filter(Chr == chr)
  completo <- agregar_gaps_telomeros(tmp, chr_sizes[chr], chr)
  datos_completo <- bind_rows(datos_completo, completo)
}

plot_data <- data.frame()

for(i in seq_along(cromosomas)) {
  chr_actual <- cromosomas[i]
  chr_datos <- datos_completo %>% filter(Chr == chr_actual)
  y_base <- (i - 1) * 7

  hap1 <- chr_datos %>%
    mutate(
      ymin = y_base,
      ymax = y_base + 1.4,
      fill_color = colorhap1
    )

  hap2 <- chr_datos %>%
    mutate(
      ymin = y_base + 2.2,
      ymax = y_base + 3.6,
      fill_color = colorhap2
    )

  plot_data <- bind_rows(plot_data, hap1, hap2)
}

grafico <- ggplot(plot_data) +
  geom_rect(
    aes(
      xmin = Inicio,
      xmax = Fin,
      ymin = ymin,
      ymax = ymax,
      fill = fill_color
    )
  ) +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(
    title = paste0(
      "Local Ancestry Inference: ",
      etiqueta_set,
      " (Individual ",
      individuo,
      ")"
    ),
    x = "Position (Mb)",
    y = ""
  ) +
  scale_y_continuous(limits = c(0, 160)) +
  scale_x_continuous(
    limits = c(-20 * 1e6, 250 * 1e6),
    breaks = seq(0 * 1e6, 250 * 1e6, by = 10 * 1e6),
    labels = seq(0, 250, by = 10)
  )

for(i in seq_along(cromosomas)) {
  chr_actual <- cromosomas[i]
  y_base <- (i - 1) * 7

  grafico <- grafico +
    annotate(
      "text",
      x = -20000000,
      y = y_base + 1.8,
      label = paste0(chr_actual, " - ", etiqueta_set),
      hjust = 0,
      size = 2.5
    )
}

grafico <- grafico +
  annotate(
    "rect",
    xmin = 190e6, xmax = 250e6, ymin = 125, ymax = 155,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  annotate("rect", xmin = 195e6, xmax = 202e6, ymin = 149, ymax = 152, fill = "#1B4F72") +
  annotate("text", x = 204e6, y = 150.5, label = "Europe", hjust = 0, size = 3.5, fontface = "italic") +
  annotate("rect", xmin = 195e6, xmax = 202e6, ymin = 143, ymax = 146, fill = "#C62828") +
  annotate("text", x = 204e6, y = 144.5, label = "Native American indigenous", hjust = 0, size = 3.5, fontface = "italic") +
  annotate("rect", xmin = 195e6, xmax = 202e6, ymin = 137, ymax = 140, fill = "#F4D03F") +
  annotate("text", x = 204e6, y = 138.5, label = "Sub-Saharian Africa", hjust = 0, size = 3.5, fontface = "italic") +
  annotate("rect", xmin = 195e6, xmax = 202e6, ymin = 131, ymax = 134, fill = "black") +
  annotate("text", x = 204e6, y = 132.5, label = "Unknown", hjust = 0, size = 3.5, fontface = "italic") +
  annotate("rect", xmin = 195e6, xmax = 202e6, ymin = 126, ymax = 129, fill = "gray") +
  annotate("text", x = 204e6, y = 127.5, label = "GAP", hjust = 0, size = 3.5, fontface = "italic")

nombre_salida <- file.path(
  base_dir,
  paste0("LocalAncestry_", set_datos, "_", individuo, ".tiff")
)

ggsave(
  nombre_salida,
  grafico,
  width = 32,
  height = 12,
  units = "cm",
  dpi = 600,
  compression = "lzw"
)

cat("Saved:", nombre_salida, "\n")
