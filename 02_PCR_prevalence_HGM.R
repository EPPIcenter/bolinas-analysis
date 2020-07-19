library(rstan)
library(bayesplot)
library(patchwork)

set.seed(134534535)

## Estimation of PCR prevalence using the hypergeometric distribution
Stan_model_HGM <- '
data {
	
	int y;
	int sample_size;
	int pop_size;
	
	real Se;
	real Sp;
	
}
parameters {
	
	real<lower=0.0, upper=1.0> prev;
	
}
transformed parameters {
	
	real a = (prev*Se + (1.0-prev)*(1.0-Sp)) * pop_size;
	real b = pop_size - a;
	
	real log_lik = lchoose(a,y) + lchoose(b, sample_size-y) - lchoose(pop_size, sample_size);
	
}
model {
	
	target += log_lik;
	
}
'

## Fit the model
fit_Stan_HGM_0.8 <- stan(
	model_code=Stan_model_HGM,
	data=list(y=0, sample_size=1282, pop_size=round(1282/0.8), Se=0.8, Sp=1.0))

fit_Stan_HGM_0.9 <- stan(
	model_code=Stan_model_HGM,
	data=list(y=0, sample_size=1282, pop_size=round(1282/0.9), Se=0.8, Sp=1.0))

fit_Stan_HGM_1.0 <- stan(
	model_code=Stan_model_HGM,
	data=list(y=0, sample_size=1282, pop_size=1282, Se=0.8, Sp=1.0))

options(scipen=9999)

print(fit_Stan_HGM_0.8, pars=c("prev"), digits=5)
print(fit_Stan_HGM_0.9, pars=c("prev"), digits=5)
print(fit_Stan_HGM_1.0, pars=c("prev"), digits=5)

## Plot
p_0.8 <- mcmc_trace(fit_Stan_HGM_0.8, pars=c("prev")) + theme_bw() + ggtitle("0.8")
p_0.9 <- mcmc_trace(fit_Stan_HGM_0.9, pars=c("prev")) + theme_bw() + ggtitle("0.9")
p_1.0 <- mcmc_trace(fit_Stan_HGM_1.0, pars=c("prev")) + theme_bw() + ggtitle("1.0")

p_0.8 + p_0.9 + p_1.0
