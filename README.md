# Validación de Brújulas Digitales vs. Brunton en Escenarios Experimentales

Este repositorio contiene el script oficial en R y las directrices de archivos para la réplica estadística del artículo científico: **"Variabilidad en la medición de orientaciones estructurales mediante brújula geológica y dispositivos móviles: análisis en escenarios mono-operador y multi-operador"**.

El algoritmo ejecuta un análisis metrológico multivariante completo de forma automatizada (remoción de valores atípicos mediante el protocolo de Tukey, cálculo de errores y pruebas de paridad distributiva de sensores MEMS).

## 📊 Estructura de Archivos y Preparación de Datos

Para ejecutar el análisis de manera independiente para cada escenario del artículo, usted debe preparar **dos archivos CSV distintos** en su directorio de trabajo, asegurando que las columnas coincidan exactamente con lo requerido por el script. Revise los achivos utilizados en el artículo.

### 1. Archivos de Entrada Requeridos
* **`Datos_estructurales_numerados.csv`**: Corresponde al **Escenario 1 (Mono-operador)**. Contiene las mediciones secuenciales controladas sobre 10 planos fijos (un gran total de 50 registros en la matriz base) para aislar y evaluar la precisión intrínseca del hardware.
* **`Datos_estructurales_aleatorios.csv`**: Corresponde al **Escenario 2 (Multi-operador)**. Contiene la base de datos con los 80 registros aleatorios capturados por múltiples profesionales en campo para evaluar la reproducibilidad e incertidumbre operativa inter-usuario.

### 2. Formato Obligatorio de las Columnas (CSV)
Ambos archivos deben tener la misma estructura matricial. El script normaliza automáticamente los encabezados a minúsculas, por lo que solo debe asegurarse de incluir estas columnas:

* `modelo.celular`: Marca y modelo comercial del dispositivo móvil evaluado (ej. *iPhone 13 Pro Max*, *Oukitel G1*, *Infinix Note 50S*, *Redmi Note 8 Pro*, *Xiaomi Poco F3*).
* `dip`: Ángulo de buzamiento medido con la aplicación digital (*FieldMove Clino*).
* `dipaz`: Dirección de buzamiento (azimut) medida con la aplicación digital.
* `brunton_dip`: Ángulo de buzamiento de referencia medido con la brújula Brunton Transit.
* `brunton_dipaz`: Dirección de buzamiento de referencia medida con la brújula Brunton Transit.

---

## 🛠️ Instrucciones de Ejecución en R

El código está completamente automatizado y genera todas las salidas por sí solo. Solo debe especificar al inicio del script qué archivo y dirección desea procesar:

### Para ejecutar el Escenario Mono-operador (Figura 4):
En la sección **1. CARGA DE DATOS** de su script, defina la variable `archivo` apuntando al set de datos numerados:
```R
setwd("SU_DIRECCION_LOCAL/Clino vs Brujula")
archivo <- "Datos_estructurales_numerados.csv"

Para procesar el entorno aleatorio inter-usuario, simplemente cambie el nombre del archivo en la misma línea:

setwd("SU_DIRECCION_LOCAL/Clino vs Brujula")
archivo <- "Datos_estructurales_aleatorios.csv"

💾 Salidas Automatizadas del Script

Tras presionar Run, el algoritmo procesará el archivo asignado y guardará automáticamente en su directorio de trabajo los siguientes productos listos para el formato editorial del artículo:

    Tabla_con_outliers.csv y Tabla_sin_outliers.csv: Reportes con métricas metrológicas clave (N, MAE, RMSE, Bias y el Coeficiente de Concordancia CCC de Lin para ambas componentes).

    Dunn_Dip_...csv y Dunn_Dir_...csv: Matrices de la prueba no paramétrica post-hoc de Dunn con ajuste de Bonferroni para evaluar diferencias distributivas significativas entre dispositivos.

    FIGURA_COMPARATIVA_TOTAL.png: Mosaico gráfico compuesto de alta resolución (18×12 pulgadas) que integra los paneles de Bland-Altman (A), Concordancia (B), Diagramas de caja de variabilidad del error (C) y Curvas de densidad de probabilidad (D).

📄 Licencia

Este software científico se distribuye bajo la licencia Creative Commons Attribution 4.0 International (CC-BY-4.0). Se permite su libre uso, modificación y distribución condicionado a la citación formal del artículo original y de este repositorio.
