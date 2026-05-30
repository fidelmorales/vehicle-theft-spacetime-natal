# Assessing the Impact of Public Policies on Vehicle Theft Incidence in Natal

This repository contains the supplementary material, source code, and reproducibility information for the manuscript:

**Assessing the Impact of Public Policies on Vehicle Theft Incidence in Natal: A Spatial-Temporal Model Using Non-Homogeneous Poisson Processes**

## Repository contents

This repository includes:

- Supplementary material for the manuscript;
- R code used to reproduce the analyses;
- Instructions for loading the data;
- Instructions for reproducing the main tables and figures;
- Links to archived MCMC outputs used in the analysis.

## Data

The data used in the manuscript are available through the R package `STPoissonSS`.

To install the package, use:

```r
install.packages("devtools")
devtools::install_github("Projeto-CNPq-Clima/STPoissonSS", force = TRUE)
library(STPoissonSS)
