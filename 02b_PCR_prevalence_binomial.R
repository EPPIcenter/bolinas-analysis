library(rstan)
library(bayesplot)
library(patchwork)

set.seed(134534535)

## Estimation of PCR prevalence using the binomial distribution
Stan_model_binomial <- '
data {
	
	int y;
	int sample_size;
	
	real Se;
	real Sp;
	
}
parameters {
	
	real<lower=0.0, upper=1.0> prev;
	
}
transformed parameters {
	
	real p_sample = prev*Se + (1.0-prev)*(1.0-Sp);
	
}
model {
	
	y ~ binomial(sample_size, p_sample);
	
}
'

## Fit the model
fit_Stan_binomial_0.8 <- stan(
	model_code=Stan_model_binomial,
	data=list(y=0, sample_size=1282, Se=0.8, Sp=1.0))

fit_Stan_binomial_0.9 <- stan(
	model_code=Stan_model_binomial,
	data=list(y=0, sample_size=1282, Se=0.8, Sp=1.0))

fit_Stan_binomial_1.0 <- stan(
	model_code=Stan_model_binomial,
	data=list(y=0, sample_size=1282, Se=0.8, Sp=1.0))

options(scipen=9999)

print(fit_Stan_binomial_0.8, pars=c("prev"), digits=5)
print(fit_Stan_binomial_0.9, pars=c("prev"), digits=5)
print(fit_Stan_binomial_1.0, pars=c("prev"), digits=5)

## Plot
p_0.8 <- mcmc_trace(fit_Stan_binomial_0.8, pars=c("prev")) + theme_bw() + ggtitle("0.8")
p_0.9 <- mcmc_trace(fit_Stan_binomial_0.9, pars=c("prev")) + theme_bw() + ggtitle("0.9")
p_1.0 <- mcmc_trace(fit_Stan_binomial_1.0, pars=c("prev")) + theme_bw() + ggtitle("1.0")

p_0.8 + p_0.9 + p_1.0
