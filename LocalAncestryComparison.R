#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

# =========================
# ARGUMENTS
# =========================

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1 || length(args) > 2) {
  stop("Usage: Rscript LocalAncestryComparison.R <individual_id> [base_directory]")
}

individuo <- args[1]
base_dir <- ifelse(length(args) == 2, args[2], "/media/daniel/Espacio/NuevosBED")

cat("Processing individual:", individuo, "\n")
cat("Base directory:", base_dir, "\n")

# =========================
# FILES
# =========================

archivo_LC <- file.path(
  base_dir,
  "LC0.9_GAP",
  paste0("InputG", individuo, ".tsv")
)

archivo_HC <- file.path(
  base_dir,
  "HC0.9_GAP",
  paste0("InputG", individuo, ".tsv")
)

# =========================
# CHECK FILES
# =========================

if(!file.exists(archivo_LC)) {
  stop(paste("File does not exist:", archivo_LC))
}

if(!file.exists(archivo_HC)) {
  stop(paste("File does not exist:", archivo_HC))
}

# =========================
# READ DATA
# =========================

datos_LC <- read.csv(archivo_LC, sep = "\t")
datos_HC <- read.csv(archivo_HC, sep = "\t")

# =========================
# CHROMOSOMES
# =========================

cromosomas <- paste0("chr", 1:22)

# =========================
# hg38 CHROMOSOME SIZES
# =========================

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

# =========================
# GAP FUNCTION
# =========================

agregar_gaps_telomeros <- function(df, chr_size) {
  
  df <- df %>%
    arrange(Inicio)
  
  resultado <- data.frame()
  
  # GAP inicial
  if(df$Inicio[1] > 0) {
    
    gap_ini <- data.frame(
      Chr = df$Chr[1],
      Inicio = 0,
      Fin = df$Inicio[1],
      colorhap1 = "gray",
      colorhap2 = "gray"
    )
    
    resultado <- bind_rows(resultado, gap_ini)
  }
  
  # recorrer segmentos
  for(i in 1:(nrow(df)-1)) {
    
    resultado <- bind_rows(resultado, df[i, ])
    
    if(df$Fin[i] < df$Inicio[i+1]) {
      
      gap <- data.frame(
        Chr = df$Chr[i],
        Inicio = df$Fin[i],
        Fin = df$Inicio[i+1],
        colorhap1 = "gray",
        colorhap2 = "gray"
      )
      
      resultado <- bind_rows(resultado, gap)
    }
  }
  
  # último segmento
  resultado <- bind_rows(resultado, df[nrow(df), ])
  
  # GAP final
  if(df$Fin[nrow(df)] < chr_size) {
    
    gap_fin <- data.frame(
      Chr = df$Chr[1],
      Inicio = df$Fin[nrow(df)],
      Fin = chr_size,
      colorhap1 = "gray",
      colorhap2 = "gray"
    )
    
    resultado <- bind_rows(resultado, gap_fin)
  }
  
  return(resultado)
}

# =========================
# COMPLETE LC
# =========================

datos_LC_completo <- data.frame()

for(chr in cromosomas) {
  
  tmp <- datos_LC %>%
    filter(Chr == chr)
  
  completo <- agregar_gaps_telomeros(
    tmp,
    chr_sizes[chr]
  )
  
  datos_LC_completo <- bind_rows(
    datos_LC_completo,
    completo
  )
}

# =========================
# COMPLETE HC
# =========================

datos_HC_completo <- data.frame()

for(chr in cromosomas) {
  
  tmp <- datos_HC %>%
    filter(Chr == chr)
  
  completo <- agregar_gaps_telomeros(
    tmp,
    chr_sizes[chr]
  )
  
  datos_HC_completo <- bind_rows(
    datos_HC_completo,
    completo
  )
}

# =========================
# PLOT DATAFRAME
# =========================

plot_data <- data.frame()

for(i in seq_along(cromosomas)) {
  
  chr_actual <- cromosomas[i]
  
  chr_LC <- datos_LC_completo %>%
    filter(Chr == chr_actual)
  
  chr_HC <- datos_HC_completo %>%
    filter(Chr == chr_actual)
  
  y_base <- (i - 1) * 14
  
  # LC
  lc_hap1 <- chr_LC %>%
    mutate(
      ymin = y_base,
      ymax = y_base + 1.4,
      fill_color = colorhap1
    )
  
  lc_hap2 <- chr_LC %>%
    mutate(
      ymin = y_base + 2.2,
      ymax = y_base + 3.6,
      fill_color = colorhap2
    )
  
  # HC
  hc_hap1 <- chr_HC %>%
    mutate(
      ymin = y_base + 5.4,
      ymax = y_base + 6.8,
      fill_color = colorhap1
    )
  
  hc_hap2 <- chr_HC %>%
    mutate(
      ymin = y_base + 7.6,
      ymax = y_base + 9,
      fill_color = colorhap2
    )
  
  plot_data <- bind_rows(
    plot_data,
    lc_hap1,
    lc_hap2,
    hc_hap1,
    hc_hap2
  )
}

# =========================
# PLOT
# =========================

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
    panel.background = element_rect(fill = "white")
  ) +
  
  labs(
    title = paste0(
      "Comparison of Local Ancestry Inference: lcWGS-Imputed vs hcWGS (Individual ",
      individuo,
      ")"
    ),
    x = "Position (Mb)",
    y = ""
  ) +
  
  scale_y_continuous(limits = c(0, 310)) +
  
  scale_x_continuous(
    limits = c(-20 * 1e6, 250 * 1e6),
    breaks = seq(0 * 1e6, 250 * 1e6, by = 10 * 1e6),
    labels = seq(0, 250, by = 10)
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold"
    )
  )

# =========================
# LABELS
# =========================

for(i in seq_along(cromosomas)) {
  
  chr_actual <- cromosomas[i]
  
  y_base <- (i - 1) * 14.025
  
  grafico <- grafico +
    
    annotate(
      "text",
      x = -20000000,
      y = y_base + 1.8,
      label = paste0(chr_actual, " - lcWGS"),
      hjust = 0,
      size = 2.5
    ) +
    
    annotate(
      "text",
      x = -20000000,
      y = y_base + 6.8,
      label = paste0(chr_actual, " - hcWGS"),
      hjust = 0,
      size = 2.5
    )
}

# =========================
# LEGEND
# =========================

# Add legend inside the plot area
grafico <- grafico +
  
  # Legend frame
  annotate(
    "rect",
    xmin = 190e6, xmax = 250e6,
    ymin = 250, ymax = 300,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  
  # Europe
  annotate(
    "rect",
    xmin = 195e6, xmax = 202e6,
    ymin = 290, ymax = 293,
    fill = "#1B4F72"
  ) +
  annotate(
    "text",
    x = 204e6, y = 291.5,
    label = "Europe",
    hjust = 0,
    size = 3.5,
    fontface = "italic"
  ) +
  
  # Native American
  annotate(
    "rect",
    xmin = 195e6, xmax = 202e6,
    ymin = 282, ymax = 285,
    fill = "#C62828"
  ) +
  annotate(
    "text",
    x = 204e6, y = 283.5,
    label = "Native American indigenous",
    hjust = 0,
    size = 3.5,
    fontface = "italic"
  ) +
  
  # Sub-Saharian Africa
  annotate(
    "rect",
    xmin = 195e6, xmax = 202e6,
    ymin = 274, ymax = 277,
    fill = "#F4D03F"
  ) +
  annotate(
    "text",
    x = 204e6, y = 275.5,
    label = "Sub-Saharian Africa",
    hjust = 0,
    size = 3.5,
    fontface = "italic"
  ) +
  
  # Unknown
  annotate(
    "rect",
    xmin = 195e6, xmax = 202e6,
    ymin = 266, ymax = 269,
    fill = "black"
  ) +
  annotate(
    "text",
    x = 204e6, y = 267.5,
    label = "Unknown",
    hjust = 0,
    size = 3.5,
    fontface = "italic"
  ) +
  
  # GAP
  annotate(
    "rect",
    xmin = 195e6, xmax = 202e6,
    ymin = 258, ymax = 261,
    fill = "gray"
  ) +
  annotate(
    "text",
    x = 204e6, y = 259.5,
    label = "GAP",
    hjust = 0,
    size = 3.5,
    fontface = "italic"
  )

# =========================
# SAVE
# =========================

nombre_salida <- file.path(
  base_dir,
  paste0("LocalAncestry_Comparison_", individuo, ".tiff")
)

ggsave(
  nombre_salida,
  grafico,
  width = 32,
  height = 16,
  units = "cm",
  dpi = 600,
  compression = "lzw"
)

cat("Saved:", nombre_salida, "\n")