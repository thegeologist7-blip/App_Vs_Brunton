# ==============================================================================
# DATA REPOSITORIES AND REPRODUCIBLE SCIENCE - METROLOGICAL VALIDATION
# ==============================================================================
# Title: R script for metrological validation and multivariate statistical analysis 
#        of structural geology data captured via mobile devices vs Brunton compasses.
# Author: Byron Bravo-Granda
# Institution: Instituto de Investigación Geológico y Energético (IIGE), Ecuador
# Project: Structural orientation measurements validation under single and multi-operator settings.
# Contact: thegeologist7@mail.com
# Version: 1.0.0
# License: Creative Commons Attribution 4.0 International (CC-BY-4.0)
# ==============================================================================

# =========================================================
# VALIDACIÓN DE BRÚJULAS DIGITALES vs BRUNTON

# =========================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(DescTools)
library(knitr)
library(FSA)
library(rstatix)
library(tidyr)
library(rlang)
# =========================================================
# 1. CARGA DE DATOS
# =========================================================

setwd("G:/Mi unidad/Articulos/GeoLatitud/Clino vs Brujula/Datos numerados")
archivo <- "Datos_estructurales_numerados.csv"

df <- tryCatch({
  read.csv(archivo, sep = ";", header = TRUE, fileEncoding = "latin1")
}, error = function(e){
  read.csv(archivo, sep = ",", header = TRUE, fileEncoding = "latin1")
})

if(ncol(df) <= 1) df <- read.csv(archivo, sep = ",", header = TRUE)

colnames(df) <- tolower(make.names(colnames(df)))

# =========================================================
# 2. PROCESAMIENTO BASE
# =========================================================

df_analisis <- df %>%
  mutate(
    diff_dip = dip - brunton_dip,
    mean_dip = (dip + brunton_dip)/2,
    
    diff_dir = ((dipaz - brunton_dipaz + 180) %% 360) - 180,
    mean_dir = atan2(
      sin(dipaz*pi/180) + sin(brunton_dipaz*pi/180),
      cos(dipaz*pi/180) + cos(brunton_dipaz*pi/180)
    ) * 180/pi,
    
    modelo = as.factor(modelo.celular)
  )

# =========================================================
# 3. DETECCIÓN DE OUTLIERS
# =========================================================

df_flag <- df_analisis %>%
  group_by(modelo) %>%
  mutate(
    # Para Buzamiento (Dip)
    q1_dip = quantile(diff_dip, 0.25),
    q3_dip = quantile(diff_dip, 0.75),
    iqr_dip = q3_dip - q1_dip,
    outlier_dip = diff_dip < (q1_dip - 1.5 * iqr_dip) | diff_dip > (q3_dip + 1.5 * iqr_dip),
    
    # Para Dirección (Dir)
    q1_dir = quantile(diff_dir, 0.25),
    q3_dir = quantile(diff_dir, 0.75),
    iqr_dir = q3_dir - q1_dir,
    outlier_dir = diff_dir < (q1_dir - 1.5 * iqr_dir) | diff_dir > (q3_dir + 1.5 * iqr_dir),
    
    outlier_total = outlier_dip | outlier_dir
  ) %>%
  ungroup()

df_con <- df_flag
df_sin <- df_flag %>% filter(!outlier_total)

# Verificación en consola para que estés seguro:
cat("Datos originales:", nrow(df_con), "\n")
cat("Datos tras quitar outliers:", nrow(df_sin), "\n")
cat("Outliers eliminados:", nrow(df_con) - nrow(df_sin), "\n")

# =========================================================
# 4. FUNCIÓN GENERAL DE ANÁLISIS
# =========================================================

analizar_dataset <- function(data, sufijo){
  
  # -------- TABLA --------
  tabla <- data %>%
    group_by(modelo) %>%
    group_modify(~{
      d <- .x
      data.frame(
        N = nrow(d),
        MAE_Dip = mean(abs(d$diff_dip)),
        RMSE_Dip = sqrt(mean(d$diff_dip^2)),
        Bias_Dip = mean(d$diff_dip),
        MAE_Dir = mean(abs(d$diff_dir)),
        RMSE_Dir = sqrt(mean(d$diff_dir^2)),
        Bias_Dir = mean(d$diff_dir),
        CCC_Dip = DescTools::CCC(d$dip, d$brunton_dip)$rho.c[1],
        CCC_Dir = DescTools::CCC(d$dipaz, d$brunton_dipaz)$rho.c[1]
      )
    }) %>%
    ungroup() %>%
    mutate(across(where(is.numeric), ~round(.x,3)))
  
  write.csv(tabla, paste0("Tabla_", sufijo, ".csv"), row.names = FALSE)
  
  # -------- ESTADÍSTICA --------
  kw_dip <- kruskal.test(diff_dip ~ modelo, data = data)
  kw_dir <- kruskal.test(diff_dir ~ modelo, data = data)
  
  dunn_dip <- dunnTest(diff_dip ~ modelo, data = data, method="bonferroni")
  dunn_dir <- dunnTest(diff_dir ~ modelo, data = data, method="bonferroni")
  
  write.csv(dunn_dip$res, paste0("Dunn_Dip_", sufijo, ".csv"), row.names = FALSE)
  write.csv(dunn_dir$res, paste0("Dunn_Dir_", sufijo, ".csv"), row.names = FALSE)
  
  
  # =========================================================
  # FUNCIÓN DE ANÁLISIS CON ETIQUETAS TÉCNICAS PROFESIONALES
  # =========================================================
  
  plot_panel <- function(data, x_mean, y_diff, x_true, y_obs, titulo){
    
    xm <- rlang::sym(x_mean); yd <- rlang::sym(y_diff)
    xt <- rlang::sym(x_true); yo <- rlang::sym(y_obs)
    
    # --- ETIQUETAS DINÁMICAS TÉCNICAS ---
    if(grepl("Dir", titulo)){
      l_m <- "Media dirección de buzamiento (°)"; l_d <- "Diferencia angular (App-Brunton) (°)"
      l_t <- "Dirección Brunton (°)"; l_o <- "Dirección App (°)"
    } else {
      l_m <- "Media de buzamiento (°)"; l_d <- "Diferencia de buzamiento (App-Brunton) (°)"
      l_t <- "Buzamiento Brunton (°)"; l_o <- "Buzamiento App (°)"
    }
    
    bias <- mean(data[[rlang::as_string(yd)]], na.rm=TRUE)
    sd_d <- sd(data[[rlang::as_string(yd)]], na.rm=TRUE)
    
    # A. Bland-Altman
    p1 <- ggplot(data, aes(x = !!xm, y = !!yd, color = modelo)) +
      geom_point(alpha=0.6) +
      geom_hline(yintercept=c(bias, bias+1.96*sd_d, bias-1.96*sd_d),
                 linetype=c("solid","dashed","dashed")) +
      theme_bw() + theme(axis.title = element_text(size = 9)) + 
      labs(title=paste("A. Bland-Altman -", titulo), x=l_m, y=l_d, color = "Modelo celular")
    
    # B. Concordancia 
    p2 <- ggplot(data, aes(x = !!xt, y = !!yo, color = modelo)) +
      geom_point(alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
      theme_bw() + 
      theme(axis.title = element_text(size = 9)) + 
      labs(title = paste("B. Concordancia -", titulo), x = l_t, y = l_o, color = "Modelo celular")
    
    # C. Boxplot con puntos 
    p3 <- ggplot(data, aes(x = modelo, y = !!yd, fill = modelo)) +
      geom_boxplot(alpha=0.5, outlier.shape = NA) +
      geom_jitter(width=0.2, alpha=0.3, size=1.0) +
      geom_hline(yintercept=0, linetype="dotted") +
      theme_bw() + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(axis.title = element_text(size = 9)) + 
      labs(title=paste("C. Variabilidad del Error -", titulo), x="Diagramas de caja", y=l_d, fill = "Modelo celular")
    
    # D. Curvas de Densidad Suavizadas 
    p4 <- ggplot(data, aes(x = !!yd, fill = modelo, color = modelo)) +
      geom_density(alpha = 0.3, linewidth = 0.8) +
      theme_bw() + 
      theme(
        axis.title = element_text(size = 9),
        legend.key = element_rect(fill = "white", color = NA) # Limpia el recuadro gris de fondo
      ) + 
      labs(title = paste("D. Distribución del Error -", titulo), x = l_d, y = "Densidad", fill = "Modelo celular") +
      # El truco maestro: forzamos a la leyenda de 'fill' a heredar el 'color' del borde y congelamos el alpha
      guides(
        color = "none",
        fill = guide_legend(override.aes = list(color = scales::hue_pal()(length(levels(data$modelo))), alpha = 0.3, linewidth = 0.8))
      )
    
    return((p1 + p2) / (p3 + p4))
  }
  # -------- PANELES --------
  
  panel_dip <- plot_panel(data, "mean_dip","diff_dip","brunton_dip","dip","Dip")
  panel_dir <- plot_panel(data, "mean_dir","diff_dir","brunton_dipaz","dipaz","DipDir")
  
  ggsave(paste0("FIG_Dip_", sufijo, ".png"), panel_dip, width=14, height=10)
  ggsave(paste0("FIG_Dir_", sufijo, ".png"), panel_dir, width=14, height=10)
  
  return(list(dip=panel_dip, dir=panel_dir, kw_dip=kw_dip$p.value, kw_dir=kw_dir$p.value))
  
} 
# =========================================================
# 5. EJECUCIÓN DOBLE
# =========================================================

res_con <- analizar_dataset(df_con, "con_outliers")
res_sin <- analizar_dataset(df_sin, "sin_outliers")

# =========================================================
# 6. PANEL COMPARATIVO FINAL CON ENCABEZADOS DE BLOQUE Y COMPONENTES
# =========================================================

# --- 1. TÍTULOS SUPERIORES (CON OUTLIERS) DIVIDIDOS Y ALINEADOS ---
tit_con_dip <- wrap_elements(panel = ggplot() + 
                               labs(title = "ANÁLISIS CON OUTLIERS: BUZAMIENTO (DIP)") + 
                               theme_void() + 
                               theme(plot.title = element_text(size = 11, face = "bold", color = "darkred", hjust = 0.5, vjust = 1)))

# Títulos de la derecha 
tit_con_dir <- wrap_elements(panel = ggplot(df_con, aes(x=1, y=1, color=modelo)) + 
                               labs(title = "ANÁLISIS CON OUTLIERS: DIRECCIÓN DE BUZAMIENTO (DIPDIR)") + 
                               scale_color_discrete(guide = "none") + # Genera el espacio de la leyenda en invisible
                               theme_void() + 
                               theme(plot.title = element_text(size = 11, face = "bold", color = "darkred", hjust = 0.5, vjust = 1)))

# Títulos superiores unidos horizontalmente
fila_tit_con <- tit_con_dip | tit_con_dir
graf_con <- res_con$dip | res_con$dir


# --- 2. TÍTULOS INTERMEDIOS (SIN OUTLIERS) DIVIDIDOS Y ALINEADOS ---
tit_sin_dip <- wrap_elements(panel = ggplot() + 
                               labs(title = "ANÁLISIS SIN OUTLIERS: BUZAMIENTO (DIP)") + 
                               theme_void() + 
                               theme(plot.title = element_text(size = 11, face = "bold", color = "darkblue", hjust = 0.5, vjust = 1)))

tit_sin_dir <- wrap_elements(panel = ggplot(df_sin, aes(x=1, y=1, color=modelo)) + 
                               labs(title = "ANÁLISIS SIN OUTLIERS: DIRECCIÓN DE BUZAMIENTO (DIPDIR)") + 
                               scale_color_discrete(guide = "none") + # Genera el espacio de la leyenda en invisible
                               theme_void() + 
                               theme(plot.title = element_text(size = 11, face = "bold", color = "darkblue", hjust = 0.5, vjust = 1)))

# Títulos intermedios unidoshorizontalmente
fila_tit_sin <- tit_sin_dip | tit_sin_dir
graf_sin <- res_sin$dip | res_sin$dir


# --- 3. COMPOSICIÓN MATRICIAL FINAL ---
# Mantenemos escala original de 18x12 y pesos de franja exactos
fig_final <- fila_tit_con / graf_con / fila_tit_sin / graf_sin + 
  plot_layout(heights = c(0.06, 1, 0.06, 1))

# Guardamos el archivo con  dimensiones s exactas
ggsave("FIGURA_COMPARATIVA_TOTAL.png", fig_final, width=18, height=12)

cat("\nANÁLISIS CON Y SIN OUTLIERS COMPLETO\n")
