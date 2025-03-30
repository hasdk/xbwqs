# R/xbwqs
R package implementing Bayesian Weighted Quantile Sum Regression (BWQS) models and their extensions. The current version of the package contains three core functions:
  1. `bwqs()`: BWQS for evaluating the association between exposure mixtures and a single outcome of interest which can be continuous (Gaussian), binary (Bernoulli), or count (Poisson).
  2. `dbwqs()`: Dirichlet BWQS for evaluating the association between exposure mixtures and a compositional (multivariate) outcome consisting of proportions summing to 1. 
  3. `hbwqs()`: Hiearchical BWQS for evaluating the association between exposure mixtures and a single outcome of interest across multiple cohorts specified by the user. The same outcome types as `bwqs()` are supported.
  
# Installation 
Install the most recent version of `xbwqs` from GitHub via the `remotes` package:
<code>
library(remotes)
remotes::install_github('hasdk/xbwqs')
</code>

# References

Please cite the references below when using the `xbwqs` package:

- Saddiki H, Warren JL, Lesseur C, Colicino E. “Compositional outcomes and environmental mixtures: the Dirichlet Bayesian Weighted Quantile Sum Regression” (2025). arXiv preprint; doi.org/10.48550/arXiv.2503.21428.
- Colicino E, Ascari R, Saddiki H, Mercedes-Nieves F, Pedretti NF, Huddleston K, Wright RO, Wright RJ. "Cross-cohort mixture analysis: a data integration approach with applciations on gestational age and DNA-methylation derived gestational age" (2024). Biometrical Journal; doi: 10.1002/bimj.202300270. PMID: 39473136.
- Colicino E, Pedretti NF, Busgang SA, Gennings C. "Per- and poly-fluoroalkyl substances and bone mineral density: Results from the Bayesian weighted quantile sum regression" (2020). Environmental Epidemiology; doi: 10.1097/EE9.0000000000000092. PMID: 32613152; PMCID: PMC7289141.
