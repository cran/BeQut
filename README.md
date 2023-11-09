  <!-- badges: start -->
  [![R-CMD-check](https://github.com/AntoineBbi/BeQut/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AntoineBbi/BeQut/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

# BeQut

**BeQut** is a R-package for Bayesian estimation of quantile regression mixed models. Based on the asymmetric Laplace distribution, it also allows to estimate joint models for longitudinal and time-to-event data, linear mixed effects models and simple linear models.

Yang, M., Luo, S., & DeSantis, S. (2019). Bayesian quantile regression joint models: Inference and dynamic predictions. Statistical Methods in Medical Research, 28(8), 2524–2537. https://doi.org/10.1177/0962280218784757

To try the current development version from github, use:

```{r} 
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")}
devtools::install_github("AntoineBbi/BeQut")
 ```
**Warning:** **BeQut** package requires JAGS software (http://mcmc-jags.sourceforge.net/). 
