library(rstan)
library(bayesplot)

set.seed(134534535)

## Estimation of sero-prevalence separately by assay, using the hypergeometric (HGM) distribution
Stan_model_separate_HGM <- '
data {
	
	int y_sample_Abbott;
	int n_sample_Abbott;
	
	int y_sample_ELISA;
	int n_sample_ELISA;
	
	int pop_size;
	
	int y_Se_Abbott;
	int n_Se_Abbott;
	
	int y_Sp_Abbott;
	int n_Sp_Abbott;
	
	int y_Se_ELISA;
	int n_Se_ELISA;
	
	int y_Sp_ELISA;
	int n_Sp_ELISA;
	
}
parameters {
	
	real<lower=0.0, upper=0.01> prev_Abbott;
	real<lower=0.0, upper=0.01> prev_ELISA;
	
	real logit_Se_Abbott;
	real logit_Sp_Abbott;
	
	real logit_Se_ELISA;
	real logit_Sp_ELISA;
	
}
transformed parameters {
	
	real<lower=0.0, upper=1.0> Se_Abbott = inv_logit(logit_Se_Abbott);
	real<lower=0.0, upper=1.0> Sp_Abbott = inv_logit(logit_Sp_Abbott);
	
	real<lower=0.0, upper=1.0> Se_ELISA = inv_logit(logit_Se_ELISA);
	real<lower=0.0, upper=1.0> Sp_ELISA = inv_logit(logit_Sp_ELISA);
	
	real a_Abbott = (prev_Abbott*Se_Abbott + (1.0-prev_Abbott)*(1.0-Sp_Abbott)) * pop_size;
	real b_Abbott = pop_size - a_Abbott;
	
	real a_ELISA = (prev_ELISA*Se_ELISA + (1.0-prev_ELISA)*(1.0-Sp_ELISA)) * pop_size;
	real b_ELISA = pop_size - a_ELISA;
	
	real log_lik_Abbott = lchoose(a_Abbott,y_sample_Abbott) + lchoose(b_Abbott, n_sample_Abbott-y_sample_Abbott) - lchoose(pop_size, n_sample_Abbott);
	
	real log_lik_ELISA = lchoose(a_ELISA,y_sample_ELISA) + lchoose(b_ELISA, n_sample_ELISA-y_sample_ELISA) - lchoose(pop_size, n_sample_ELISA);
	
}
model {
	
	// Sero-prevalence
	target += log_lik_Abbott;
	target += log_lik_ELISA;
	
	// Test performance characteristics
	target += binomial_lpmf(y_Se_Abbott | n_Se_Abbott, Se_Abbott);
	target += binomial_lpmf(y_Sp_Abbott | n_Sp_Abbott, Sp_Abbott);
	
	target += binomial_lpmf(y_Se_ELISA | n_Se_ELISA, Se_ELISA);
	target += binomial_lpmf(y_Sp_ELISA | n_Sp_ELISA, Sp_ELISA);
	
	// Priors
	target += normal_lpdf(logit_Se_Abbott | 5, 2);
	target += normal_lpdf(logit_Sp_Abbott | 5, 2);
	target += normal_lpdf(logit_Se_ELISA | 5, 2);
	target += normal_lpdf(logit_Sp_ELISA | 5, 2);
	
}
generated quantities {
	
	real<lower=0.0, upper=1.0> PPV_Abbott;
	real<lower=0.0, upper=1.0> PPV_ELISA;
	
	real<lower=0.0, upper=1.0> NPV_Abbott;
	real<lower=0.0, upper=1.0> NPV_ELISA;
	
	PPV_Abbott = (Se_Abbott*prev_Abbott)/((Se_Abbott*prev_Abbott) + (1.0-Sp_Abbott)*(1.0-prev_Abbott));
	
	PPV_ELISA = (Se_ELISA*prev_ELISA)/((Se_ELISA*prev_ELISA) + (1.0-Sp_ELISA)*(1.0-prev_ELISA));
	
	NPV_Abbott = (Sp_Abbott*(1.0-prev_Abbott))/((Sp_Abbott*(1.0-prev_Abbott)) + (1.0-Se_Abbott)*prev_Abbott);
	
	NPV_ELISA = (Sp_ELISA*(1.0-prev_ELISA))/((Sp_ELISA*(1.0-prev_ELISA)) + (1.0-Se_ELISA)*prev_ELISA);
	
}	
'

## Fit the model
fit_Stan_separate_HGM <- stan(
	model_code=Stan_model_separate_HGM,
	data=list(
		y_sample_Abbott=6, n_sample_Abbott=1255,
		y_sample_ELISA=3, n_sample_ELISA=1307,
		pop_size=as.integer(round(1307/0.8)),
		y_Se_Abbott=88, n_Se_Abbott=88,
		y_Sp_Abbott=1066, n_Sp_Abbott=1070,
		y_Se_ELISA=42, n_Se_ELISA=44,
		y_Sp_ELISA=95, n_Sp_ELISA=95),
	init=function(){return(list(prev_Abbott=0.001, prev_ELISA=0.001, logit_Se_Abbott=plogis(0.99), logit_Sp_Abbott=plogis(0.99), logit_Se_ELISA=plogis(0.99), logit_Sp_ELISA=plogis(0.99)))},
	chains=4, seed=234, iter=200000, thin=100)

print(fit_Stan_separate_HGM, digits=4)

color_scheme_set("brewer-Spectral")
mcmc_trace(fit_Stan_separate_HGM) + theme_bw()

color_scheme_set("brightblue")
mcmc_pairs(fit_Stan_separate_HGM, pars=c("prev_Abbott", "prev_ELISA", "Se_Abbott", "Se_ELISA", "Sp_Abbott", "Sp_ELISA", "a_Abbott", "b_Abbott", "a_ELISA", "b_ELISA"))

save(fit_Stan_separate_HGM, file="03.RData")
