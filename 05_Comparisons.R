library(rstan)
library(bayesplot)

## Compare the models
load("03.RData")
load("03b.RData")
load("04.RData")
load("04b.RData")

## Separate fitting
windows();mcmc_hist(fit_Stan_separate_HGM) + ggtitle("Separate, HGM") + theme_bw()
windows();mcmc_hist(fit_Stan_separate_binomial) + ggtitle("Separate, binomial") + theme_bw()

windows();mcmc_trace(fit_Stan_separate_HGM) + ggtitle("Separate, HGM") + theme_bw()
windows();mcmc_trace(fit_Stan_separate_binomial) + ggtitle("Separate, binomial") + theme_bw()

## Joint fitting
windows();mcmc_hist(fit_Stan_joint_multivariate_HGM) + ggtitle("Joint, multivariate HGM") + theme_bw()
windows();mcmc_hist(fit_Stan_joint_multinomial) + ggtitle("Joint, multinomial") + theme_bw()

windows();mcmc_trace(fit_Stan_joint_multivariate_HGM) + ggtitle("Joint, multivariate HGM") + theme_bw()
windows();mcmc_trace(fit_Stan_joint_multinomial) + ggtitle("Joint, multinomial") + theme_bw()
