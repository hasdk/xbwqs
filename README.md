# R/xbwqs
R package implementing Bayesian Weighted Quantile Sum Regression (BWQS) models and their variants. The current version of the package contains three core functions:
  1. `bwqs()': BWQS for evaluating the association between exposure mixtures and a single outcome of interest which can be continuous (Gaussian), binary (bernoulli), or count (Poisson).
  2. 'dbwqs()': Dirichlet BWQS for evaluating the association exposure mixtures and a compositional (multivariate) outcome consisting of proportions summing to 1. 
  3. 'hbwqs()': Hiearchical BWQS for evaluating the association between between exposure mixtures and a single outcome of interest across multiple cohorts specified by the user. The same outcome types as `bwqs()` are supported.
  
# Installation 

# DBWQS 
