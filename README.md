# rEGA <img src="man/figures/logo.png" align="right" height="120" alt="" />

[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Open Science: OSF Ready](https://img.shields.io/badge/Open%20Science-SHA256%20Audited-orange.svg)](https://osf.io/)

Paquete de R (`rEGA`) ombojoajúva red psicométrica, operador resolvente Neumann-Markov, descomposición exacta de varianza nodal (`h2 + beta2 + psi2 = 1` / `100%`) ha coeficiente de fiabilidad estructural (`Gamma_c`, `Gamma_t`, `Psi_G`, `SEM`).

---

## 1. Instalación (`GitHub`)

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Huerta-Huacchahuilca/rEGA")
```

---

## 2. Estimación del Modelo (`rEGA`, `seed = 11`)

```r
library(rEGA)

set.seed(11)
Sigma <- matrix(0.28, 6, 6)
Sigma[1:3, 1:3] <- 0.68
Sigma[4:6, 4:6] <- 0.68
diag(Sigma) <- 1.0
dat <- as.data.frame(MASS::mvrnorm(n = 300, mu = rep(0, 6), Sigma = Sigma))
colnames(dat) <- paste0("Item_", 1:6)

fit <- rEGA(
  data       = dat,
  items      = 1:6,
  cor.method = "cor_auto",
  algorithm  = "walktrap",
  plot.rEGA  = TRUE,
  seed       = 11
)

print(fit)
summary(fit)
```

---

## 3. Visualización (`plot.rEGA`: `both`, `separate`, `polychrome`, `colorblind`, `fr`, `circle`)

```r
# Figura combinada (Red + Tríada Nodal en % con umbral mínimo de etiqueta del 10%)
plot(fit, type = "both", palette = "polychrome", layout = "fr", show.values = TRUE, label.min.pct = 10)

# Figuras separadas en paneles independientes
plot(fit, type = "separate", palette = "polychrome")

# Configuración para daltónicos (Okabe-Ito) y disposición circular
plot(fit, type = "both", palette = "colorblind", layout = "circle")

# Gráficos individuales
plot(fit, type = "network", palette = "polychrome", layout = "fr")
plot(fit, type = "variance", palette = "polychrome", show.values = TRUE, label.min.pct = 10)
```

---

## 4. Inferencia Bootstrap con Barra de Progreso y ETA (`boot.rEGA`)

```r
boot_fit <- boot.rEGA(
  rEGA.obj       = fit,
  iter           = 500,
  type           = "parametric",
  ncores         = 2,
  seed           = 11,
  progress       = TRUE,
  plot.boot.rEGA = TRUE
)

print(boot_fit)
summary(boot_fit)

# Opciones de visualización para boot.rEGA
plot(boot_fit, type = "both", palette = "polychrome")
plot(boot_fit, type = "separate", palette = "colorblind")
plot(boot_fit, type = "all", palette = "polychrome")
plot(boot_fit, type = "stability", palette = "colorblind")
```

---

## 5. Potencia Estadística y Auditoría Criptográfica (`rEGA.power` & `rEGA.audit`)

```r
# Planificación de tamaño muestral con corrección por atenuación
rEGA.power(rEGA.obj = fit, p.target = c(6, 9, 12), r.true = 0.30)

# Exportación de libro Excel de auditoría con hash SHA-256
rEGA.audit(rEGA.obj = fit, boot.obj = boot_fit, file = "rEGA_Audit_Report.xlsx")
```

---

## Autor y Contacto

* **César Huerta-Huacchahuilca**  
  *Correo:* [andrehuerta1411@gmail.com](mailto:andrehuerta1411@gmail.com)  
  *GitHub:* [@Huerta-Huacchahuilca](https://github.com/Huerta-Huacchahuilca)  
  *ORCID:* [0009-0000-8906-2678](https://orcid.org/0009-0000-8906-2678)  
  *Issues:* [https://github.com/Huerta-Huacchahuilca/rEGA/issues](https://github.com/Huerta-Huacchahuilca/rEGA/issues)

---

## Citación (APA 7 & BibTeX)

> Huerta-Huacchahuilca, C. (2026). *rEGA: Reliability Exploratory Graph Analysis* (R package version 0.1.0) [Computer software]. GitHub. https://github.com/Huerta-Huacchahuilca/rEGA

```bibtex
@manual{huerta2026rega,
  title  = {rEGA: Reliability Exploratory Graph Analysis},
  author = {Huerta-Huacchahuilca, C{\'e}sar},
  year   = {2026},
  note   = {R package version 0.1.0},
  url    = {https://github.com/Huerta-Huacchahuilca/rEGA}
}
```
