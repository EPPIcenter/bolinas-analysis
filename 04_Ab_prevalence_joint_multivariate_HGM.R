library(rstan)
library(bayesplot)

set.seed(134453568)

## Estimation of seroprevalence assuming conditional independence between assay, using the multivariate hypergeomtric distribution
Stan_model_joint_multivariate_HGM <- '
data {
	
	// This is ordered as (Abbott, ELISA): (+,+), (+,-), (-,+), (-,-)
	int y_sample[4];
	
	// This is ordered as (Abbott, ELISA): (-,NA), (NA,-)
	int y_sample_singletons[2];
	
	// For hypergeometric
	int sample_size_tested_both;
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
	
	real<lower=0.0, upper=0.01> prev;
	
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
	
	simplex[4] p_sample;
	
	p_sample[1] = prev*Se_Abbott*Se_ELISA + (1.0-prev)*(1.0-Sp_Abbott)*(1.0-Sp_ELISA);
	p_sample[2] = prev*Se_Abbott*(1.0-Se_ELISA) + (1.0-prev)*(1.0-Sp_Abbott)*Sp_ELISA;
	p_sample[3] = prev*(1.0-Se_Abbott)*Se_ELISA + (1.0-prev)*Sp_Abbott*(1.0-Sp_ELISA);
	p_sample[4] = prev*(1.0-Se_Abbott)*(1.0-Se_ELISA) + (1.0-prev)*Sp_Abbott*Sp_ELISA;
	
}
model {
	
	// Sero-prevalence contributions from samples tested on both
	target += lchoose(p_sample[1] * pop_size, y_sample[1]) +
	lchoose(p_sample[2] * pop_size, y_sample[2]) +
	lchoose(p_sample[3] * pop_size, y_sample[3]) +
	lchoose(p_sample[4] * pop_size, y_sample[4]) - lchoose(pop_size, sample_size_tested_both);
	
	// Sero-prevalence contributions from singletons
	target += binomial_lpmf(0 | y_sample_singletons[1], prev*Se_Abbott + (1.0-prev)*(1.0-Sp_Abbott));
	target += binomial_lpmf(0 | y_sample_singletons[2], prev*Se_ELISA + (1.0-prev)*(1.0-Sp_ELISA));
	
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
	
	real PPV_Abbott_pos_ELISA_pos;
	real PPV_Abbott_pos_ELISA_neg;
	real PPV_Abbott_neg_ELISA_pos;
	
	PPV_Abbott_pos_ELISA_pos = (Se_Abbott*Se_ELISA*prev) / ((Se_Abbott*Se_ELISA*prev) + ((1.0-Sp_Abbott)*(1.0-Sp_ELISA)*(1.0-prev)));
	
	PPV_Abbott_pos_ELISA_neg = (Se_Abbott*(1.0-Se_ELISA)*prev) / ((Se_Abbott*(1.0-Se_ELISA)*prev) + ((1.0-Sp_Abbott)*Sp_ELISA*(1.0-prev)));
	
	PPV_Abbott_neg_ELISA_pos = ((1.0-Se_Abbott)*Se_ELISA*prev) / (((1.0-Se_Abbott)*Se_ELISA*prev) + (Sp_Abbott*(1.0-Sp_ELISA)*(1.0-prev)));
	
}
'

## Fit the model, with updated data
fit_Stan_joint_multivariate_HGM <- stan(
	model_code=Stan_model_joint_multivariate_HGM,
	data=list(
		## This is ordered as (Abbott, ELISA): (+,+), (+,-), (-,+), (-,-)
		y_sample=c(1,5,2,1247),
		## This is ordered as (Abbott, ELISA): (-,NA), (NA,-)
		y_sample_singletons=c(0,52),
		## For HGM
		sample_size_tested_both=as.integer(sum(c(1,5,2,1247))), pop_size=as.integer(round(1307/0.8)),
		## Test performance characteristics
		y_Se_Abbott=88, n_Se_Abbott=88,
		y_Sp_Abbott=1066, n_Sp_Abbott=1070,
		y_Se_ELISA=42, n_Se_ELISA=44,
		y_Sp_ELISA=95, n_Sp_ELISA=95),
	init=function(){return(list(prev=0.001, logit_Se_Abbott=plogis(0.99), logit_Sp_Abbott=plogis(0.99), logit_Se_ELISA=plogis(0.99), logit_Sp_ELISA=plogis(0.99)))},
	chains=4, seed=3453, iter=200000, thin=100)

print(fit_Stan_joint_multivariate_HGM, digits=4)

color_scheme_set("brewer-Spectral")
mcmc_trace(fit_Stan_joint_multivariate_HGM) + theme_bw()

color_scheme_set("brightblue")
mcmc_pairs(fit_Stan_joint_multivariate_HGM, pars=c("prev", "logit_Se_Abbott", "logit_Se_ELISA", "logit_Sp_Abbott", "logit_Sp_ELISA", "p_sample[1]", "p_sample[2]", "p_sample[3]", "p_sample[4]"))

save(fit_Stan_joint_multivariate_HGM, file="04.RData")
