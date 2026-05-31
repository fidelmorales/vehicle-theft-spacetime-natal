# README_REPRODUCIBILITY

## Modeling Vehicle-Theft Dynamics Across Administrative Periods in Natal, Brazil

This document provides instructions for reproducing the analyses reported in the manuscript.

---

# Step 1 – Install the STPoissonSS package

```r
install.packages("devtools")
devtools::install_github("Projeto-CNPq-Clima/STPoissonSS")
library(STPoissonSS)
```

---

# Step 2 – Download Archived MCMC Outputs

Download the archived MCMC outputs from:

https://doi.org/10.5281/zenodo.17011948

Place the files in the working directory.

---

# Step 3 – Run the Reproducibility Script

Open and execute:

```r
source("reproducibility_code.R")
```

The script contains the code used to reproduce:

* model estimation;
* convergence diagnostics;
* posterior summaries;
* Table 1;
* Table 2;
* Figure 2;
* Figure 3;
* Figure 4;
* supplementary analyses.

---

# Additional Information

Detailed descriptions of the dataset, model specification, prior distributions, MCMC implementation, and convergence assessment are provided in:

```text
SupplementaryMaterial.pdf
```
