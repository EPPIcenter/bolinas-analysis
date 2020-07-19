library(ggplot2)
library(dplyr)

n = 1282 ## Sample size
x = 0    ## Number of observed cases

## Testing the upper bound on the true number of infections

## Code adapted from: https://www.nature.com/articles/s41467-018-06657-5 
PFFI_ST_KCount <- function(K,Se,Sp,N,n,x,alpha,beta){
	
	# probability of null hypothesis
	d = K # diseased in population
	# d = floor(p*N) # diseased in population - safer as above!
	P_null = 0
	for (x1 in 0:x){ # summation of P(T+=x)
	c = 0
	for (y in 0:d){ # outer summation of hypergeometric
	a = dhyper(x=y, m=d, n=N-d, k=n) ## ST: use the distribution
	# a = (choose(d,y)*choose(N-d,n-y))/choose(N,n)
	b = 0
	for (j in 0:min(x1,y)){ # inner summation of hypergeometric
	b = b + dbinom(x=j, size=y, prob=Se)*dbinom(x=x1-j, size=n-y, prob=1-Sp) ## ST: use the distribution
	# b = b + choose(y,j)*Se^j*(1-Se)^(y-j)*choose(n-y,x1-j)*(1-Sp)^(x1-j)*Sp^(n-x1-y+j)
	}
	c = c + a*b
	}
	P_null = P_null+c
	}

	# probability of alternative hypothesis
	d=0 # disease free population
	P_alt = 0
	for (x2 in x:n){ # summation of P(T+=x)
	c = 0
	for (y in 0:d){ # outer summation of hypergeometric
	a = dhyper(x=y, m=d, n=N-d, k=n) ## ST: use the distribution
	# a = (choose(d,y)*choose(N-d,n-y))/choose(N,n)
	b = 0
	for (j in 0:min(x2,y)){ # inner summation of hypergeometric
	b = b + dbinom(x=j, size=y, prob=Se)*dbinom(x=x2-j, size=n-y, prob=1-Sp) ## ST: use the distribution
	# b = b + choose(y,j)*Se^j*(1-Se)^(y-j)*choose(n-y,x2-j)*(1-Sp)^(x2-j)*Sp^(n-x2-y+j)
	}
	c = c + a*b
	}
	P_alt = P_alt+c
	}
	
	# probability of freedom = 1-P_null
	P_free = (1-P_null)
	
	# draw a conclusion regarding whether the survey indicates the population is free from infection
	if (P_null>alpha && P_alt>=(1-beta)){ ## ST
	# if (P_null>alpha && P_alt>alpha){
	conc = "insufficient evidence, sample size too small"
	}
	
	else if (P_null<=alpha && P_alt>=(1-beta)){ ## ST
	# else if (P_null<alpha && P_alt>(1-beta)){
	conc = "free from infection at the determined threshold"
	}
	
	else if (P_alt<(1-beta)){
	conc = "not free from infection"
	}
	
	# return list of probabilities ## ST: added colon for visibility
	vars = c('P_null:','P_alt:','P_free:','Decision:')
	Prob = c(P_null,P_alt,P_free,conc)
	out = data.frame(Prob,row.names = vars)
	return (out)
	
	#out2 = data.frame(value=Prob[1], row.names=vars[1])
	return(out2)
}

## Sensitivity sweep
Se_sweep <- 0.8

## Specificity sweep
Sp_sweep <- 1

## Population size sweep
N_sweep <- round(n/c(1, 0.9, 0.8))

## True prevalence sweep
K_sweep <- c(1:5)

## Set up a matrix
empty_col <- rep(NA, times=length(Se_sweep) * length(Sp_sweep) * length(N_sweep) * length(K_sweep))

dat_sweep <- data.frame(Se=empty_col, Sp=empty_col, N=empty_col, K=empty_col, p_null=empty_col, p_alt=empty_col, p_free=empty_col, decision=empty_col)

## Fix alpha and beta
dat_sweep$alpha <- 0.05
dat_sweep$beta <- 0.05

## Loop over
counter <- 1
for(i in 1:length(Se_sweep)){

for(j in 1:length(Sp_sweep)){

for(k in 1:length(N_sweep)) {

for(l in 1:length(K_sweep)) {

	## Population matrix
	dat_sweep$Se[counter] <- Se_sweep[i]
	dat_sweep$Sp[counter] <- Sp_sweep[j]
	dat_sweep$N[counter] <- N_sweep[k]
	dat_sweep$K[counter] <- K_sweep[l]

	out = PFFI_ST_KCount(
		K=dat_sweep$K[counter],
		Se=dat_sweep$Se[counter],
		Sp=dat_sweep$Sp[counter],
		N=dat_sweep$N[counter],
		n=n,
		x=x,
		alpha=dat_sweep$alpha[1],
		beta=dat_sweep$beta[1])
	
	out2 <- out %>% t() %>% as.data.frame()

	dat_sweep$p_null[counter] <- out2 %>% select("P_null:") %>% pull() %>% as.character() %>% as.numeric()
	dat_sweep$p_alt[counter] <- out2 %>% select("P_alt:") %>% pull() %>% as.character() %>% as.numeric()
	dat_sweep$p_free[counter] <- out2 %>% select("P_free:") %>% pull() %>% as.character() %>% as.numeric()
	dat_sweep$decision[counter] <- out2 %>% select("Decision:") %>% pull() %>% as.character()

	counter <- counter + 1
	rm(out, out2)

}}}}

## Plot
dat_sweep %>%
	mutate(N=factor(N, labels=c("100%", "90%", "80%"))) %>%
	ggplot(aes(x=K, y=p_null, group=N, colour=N)) +
		geom_line(size=1) +
		geom_point(size=3) +
		theme_bw(base_size=15) +
		xlab("True number of cases") +
		ylab(paste0("Pr(observing ", x, " cases | true number of cases)")) +
		theme(legend.position = c(0.8,0.85)) +
		guides(colour=guide_legend(title=paste0("% Pop. sampled\n(n=", n, ")")))
