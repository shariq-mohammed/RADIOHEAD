# RADIOHEAD

RADIOHEAD is a package to fit the model proposed in the manuscript: 
[S Mohammed](shariq-mohammed.github.io), K Bharath, S Kurtek, A Rao, [V Baladandayuthapani](bayesrx.github.io), *2020*, _RADIOHEAD: Radiogenomic Analysis Incorporating Tumor Heterogeneity in Imaging Through Densities_.

```
#####################################################
### This  code only shows an example execution of ###
### RADIOHEAD pipeline using RADIOHEAD package.   ###
#####################################################

# Change the directory to where the package is saved
setwd('C:/Users/xyz/R Package')

# install the package
devtools::install('./RADIOHEAD')

# load the package
require(RADIOHEAD)

###############
### Example ###
###############
# Note that the data is available within the package.

# Choose the pathway score for the 'EXOCYTOSIS' pathway.
y = scale(c_pathway_scores[,"EXOCYTOSIS"], center = T, scale = F)

# Choose the principal component scores based on radiomic images.
X = scale(pc_scores, center = T, scale = T)

# Choose the group membership labels for the principal components.
groups = pc_groups

# Execute the Gibbs sampling
res = groupSS(y, X, groups, Nmcmc=2e3)

# Identify the significant radiogenomic associaitons between the pathway and ...
# ... the radiomic features after a Bayesian FDR-based variable selection.
# 'associations' outputs the effect size of the significant associations.
associations = fdr_var_selection(res$b, res$x_cnames)
associations
```
