# project 2

rm(list=ls())

##### set up data ##############################################################
setwd("C:/Users/hw1220/Desktop/project")
mag <- read.csv("magnesium.csv")

# deleting unnecessary column
mag$X <- NULL

# rearrange to fit table 2 in Higgins paper
colnames(mag) <- c("trial", "name", "year", "nm", "rm", "nc", "rc")
mag <- mag[c("trial", "name", "year", "rm", "nm", "rc", "nc")]
mag[13, 1] <- 8
mag[8:12, 1] <- c(9:13)
mag <- mag[order(mag$trial),]

# examine the studies
mag$total <- mag$nm + mag$nc

sum(mag[1:14, 8]) # 4388 patients in other 14 studies
58050/4388 # how much bigger ISIS-4 is than the others

##### frequentist analysis #####################################################
# install.packages("metafor")
# install.packages("rmeta")
# install.packages("formattable")

library(metafor)
library(rmeta)
library(formattable)

### reproduce table 2

# OR in Peto FE for first 8 trials
or.p8 <- rma.peto(ai=rm, n1i=nm, ci=rc, n2i=nc,
                  data=mag[1:8,], slab=name,
                  add=1/2, to="only0", drop00=TRUE,
                  level=95, digits=2, verbose=FALSE)
summary(or.p8) # OR 0.65, CI(0.51, 0.82)

# OR in Peto FE for first 14 trials
or.p14 <- rma.peto(ai=rm, n1i=nm, ci=rc, n2i=nc,
                   data=mag[1:14,], slab=name,
                   add=1/2, to="only0", drop00=TRUE,
                   level=95, digits=2, verbose=FALSE)
summary(or.p14) # OR 0.57, CI(0.46, 0.71)

# OR in Peto FE for all trials incl ISIS-4
or.p15 <- rma.peto(ai=rm, n1i=nm, ci=rc, n2i=nc,
                   data=mag[1:15,], slab=name,
                   add=1/2, to="only0", drop00=TRUE,
                   level=95, digits=2, verbose=FALSE)
summary(or.p15) # OR 1.01, CI(0.95, 1.07)

# OR in DSL RE for first 8 trials
or.dl8 <- meta.DSL(ntrt=nm, nctrl=nc, ptrt=rm, pctrl=rc,
                   conf.level=0.95, names=name, 
                   data=mag[1:8,],  na.action=na.fail,statistic="OR")
summary(or.dl8) # OR 0.55, CI(0.34, 0.89)

# OR in DSL RE for first 14 trials
or.dl14 <- meta.DSL(ntrt=nm, nctrl=nc, ptrt=rm, pctrl=rc,
                    conf.level=0.95, names=name, 
                    data=mag[1:14,],  na.action=na.fail,statistic="OR")
summary(or.dl14) # OR 0.47, CI(0.32, 0.68)

# OR in DSL RE for all trials incl ISIS-4
or.dl15 <- meta.DSL(ntrt=nm, nctrl=nc, ptrt=rm, pctrl=rc,
                    conf.level=0.95, names=name, 
                    data=mag[1:15,],  na.action=na.fail,statistic="OR")
summary(or.dl15) # OR 0.53, CI(0.36, 0.77)

# forest plots for all 3 FE models
forest(or.p8, main="Odds Ratio for the First 8 Trials")
forest(or.p14, main="Odds Ratio for the First 14 Trials")
forest(or.p15, main="Odds Ratio for All Trials")

##### bayesian #################################################################
# install.packages("rstan", repos = "http://cran.rstudio.com", dependencies = TRUE)
library(rstan)

k <- length(mag$trial)
nc <- mag$nc
nm <- mag$nm
rc <- mag$rc
rm <- mag$rm

# Bayesian model with reference prior
model_ref <- "
data{
    int<lower=0> k;
    int<lower=0> nc[k];
    int<lower=0> nm[k];
    int<lower=0> rc[k];
    int<lower=0> rm[k];
}

parameters {
    vector[k] delta;
    real<lower=0, upper=1> pc[k];
    real<lower=0> sigma;
    real deltanew;
    real mu;
}

transformed parameters{
    real<lower=0, upper=1> pm[k];
    for (i in 1:k) {
        pm[i] <- exp(log(pc[i]/(1-pc[i]))+delta[i])/(1+exp(log(pc[i]/(1-pc[i]))+delta[i]));
    }
}

model{
    for (i in 1:k) {
        rc[i] ~ binomial(nc[i], pc[i]);
        rm[i] ~ binomial(nm[i], pm[i]);
        delta[i] ~ normal(mu, sigma);
        pc[i] ~ uniform(0, 1);
    }

    deltanew ~ normal(mu, sigma);
    mu ~ normal(0, 100);
    sigma ~ uniform(0, 100);
}"

# Bayesian model with skeptical prior
model_skp <- "
data{
    int<lower=0> k;
    int<lower=0> nc[k];
    int<lower=0> nm[k];
    int<lower=0> rc[k];
    int<lower=0> rm[k];
}

    parameters {
    vector[k] delta;
    real<lower=0, upper=1> pc[k];
    real<lower=0> sigma;
    real deltanew;
    real mu;
}

transformed parameters{
    real<lower=0, upper=1> pm[k];
    for (i in 1:k) {
        pm[i] <- exp(log(pc[i]/(1-pc[i]))+delta[i])/(1+exp(log(pc[i]/(1-pc[i]))+delta[i]));
    }
}

model{
    for (i in 1:k) {
        rc[i] ~ binomial(nc[i], pc[i]);
        rm[i] ~ binomial(nm[i], pm[i]);
        delta[i] ~ normal(mu, sigma);
        pc[i] ~ uniform(0, 1);
    }

    deltanew ~ normal(mu, sigma);
    mu ~ normal(0, 1/sqrt(32.69));
    sigma ~ uniform(0, 100);
}"

# fit reference model in stan
mod_fit_ref <- stan(model_code=model_ref, data=c("k", "nc", "nm", "rc", "rm"),
                    pars=c("delta", "mu","sigma","deltanew"),
                    iter=500000, chains=3, warmup=500, verbose=FALSE)
print(mod_fit_ref)

# fit skeptical model in stan
mod_fit_skp <- stan(model_code=model_skp, data=c("k", "nc", "nm", "rc", "rm"),
                    pars=c("delta", "mu","sigma","deltanew"), 
                    iter=500000, chains=3, warmup=500, verbose=FALSE)
print(mod_fit_skp)

# trace plots for deltanew
traceplot(mod_fit_ref, pars= "deltanew")
# plot(mod_fit_ref)

traceplot(mod_fit_skp, pars= "deltanew")
# plot(mod_fit_skp)

# histogram for reference posterior distribution
mag.sim.ref <- extract(mod_fit_ref, permuted=TRUE)
# hist(mag.sim.ref$deltanew, breaks=200, xlim=range(-4, 2),
     # xlab="Posterior distribution of log of odds ratio",
     # main="Reference prior histogram")

hist(exp(mag.sim.ref$deltanew), breaks=20000, xlim = range(0, 2), 
     xlab="Posterior distribution of odds ratio",
     main="Reference prior histogram")

# histogram for skeptical posterior distribution
mag.sim.skp <- extract(mod_fit_skp, permuted=TRUE)
# hist(mag.sim.skp$deltanew, breaks=200, xlim=range(-3, 3),
     # xlab="Posterior distribution of log of odds ratio",
     # main="Skeptical prior histogram")

hist(exp(mag.sim.skp$deltanew), breaks=200000, xlim = range(0, 4), 
     xlab="Posterior distribution of odds ratio",
     main="Skeptical prior histogram")

# statistical superiority
sum(exp(mag.sim.ref$deltanew) < 1) / length(mag.sim.ref$deltanew) # 0.90
sum(exp(mag.sim.skp$deltanew) < 1) / length(mag.sim.skp$deltanew) # 0.68

# clinical superiority
sum(exp(mag.sim.ref$deltanew) < 0.9) / length(mag.sim.ref$deltanew) # 0.87
sum(exp(mag.sim.skp$deltanew) < 0.9) / length(mag.sim.skp$deltanew) # 0.62

##### testing rstan ############################################################
set_cppo('fast')
set_cppo('debug')

install.packages("inline")
install.packages("Rcpp")
library(inline) 
library(Rcpp)
src <- ' 
std::vector<std::string> s; 
s.push_back("hello");
s.push_back("world");
return Rcpp::wrap(s);
'
hellofun <- cxxfunction(body = src, includes = '', plugin = 'Rcpp', verbose = FALSE)
cat(hellofun(), '\n') 
