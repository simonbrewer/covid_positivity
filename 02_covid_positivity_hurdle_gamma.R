## ----setup, include=FALSE-------------------------------------------------------------------

## Following BYM example @
##https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fece3.3081&file=ece33081-sup-0002-AppendixS2.txt

## -------------------------------------------------------------------------------------------
library(sf)
#library(tmap)
library(dplyr)
library(ggpubr)
library(lubridate)
library(spdep)
library(INLA)
library(INLAutils)
library(ggregplot)

inla.setOption(pardiso.license = "~/pardiso.lic")

## -------------------------------------------------------------------------------------------
load("dat2.RData")
dim(dat2)

## Time index
dat2$date1 <- as.numeric(ymd(dat2$date) - min(ymd(dat2$date)) + 1)

## SF object for plotting
dat_sf <- dat2 %>%
  filter(date1 == 1)

## -------------------------------------------------------------------------------------------
g <- inla.read.graph(filename = "./map.adj")

## Subset for testsing
# dat2 <- dat2 %>%
#   filter(date1 < 100)

## Number of samples
n2 <- nrow(dat2)
nothing1 <- nothing2 <- rep(NA, n2)

## Outcome variables
z <- as.vector(dat2$TESTSNEW7POS_CDC > 0)
y <- ifelse(z == 1, dat2$TESTSNEW7POS_CDC, NA)

#Generating 3 vectors of response variables for the joint model 
zNA = as.vector(c(z, nothing2))
yNA = as.vector(c(nothing1, y))

#Combine them in a matrix
outcome.matrix<-matrix(c(zNA, yNA), ncol=2)

#A factor vector with 2 levels indicating 2 outcome variables
mu = as.factor(c(rep(1,length(z)), rep(2,length(y))))

# Index vectors for the time effects
i.date1 <- as.integer(c(dat2$date1, nothing2))
i.date2 <- as.integer(c(nothing1, dat2$date1))

#Index vectors for the spatial effects 
i.spat1 <- c(dat2$struct, nothing2) # Binomial
i.spat2 <- c(nothing1, dat2$struct) # Rate

#The covariates 
Pop_m <- scale(dat2$Pop_m)
Pop_m1 <- c(Pop_m, nothing2) # Binomial male popn
Pop_m2 <- c(nothing1, Pop_m) # Gamma male popn

gini <- scale(dat2$gini)
gini1 <- c(gini, nothing2) # Binomial male popn
gini2 <- c(nothing1, gini) # Gamma male popn

Uninsured <- scale(dat2$Uninsured)
Uninsured1 <- c(Uninsured, nothing2) # Binomial male popn
Uninsured2 <- c(nothing1, Uninsured) # Gamma male popn

Production <- scale(dat2$Production)
Production1 <- c(Production, nothing2) # Binomial male popn
Production2 <- c(nothing1, Production) # Gamma male popn

RUCC <- dat2$RUCC
RUCC1 <- c(RUCC, nothing2) # Binomial male popn
RUCC2 <- c(nothing1, RUCC) # Gamma male popn

PCP <- scale(dat2$PCP)
PCP1 <- c(PCP, nothing2) # Binomial male popn
PCP2 <- c(nothing1, PCP) # Gamma male popn

RPL_THEMES <-  scale(dat2$RPL_THEMES)
RPL_THEMES1 <- c(RPL_THEMES, nothing2) # Binomial male popn
RPL_THEMES2 <- c(nothing1, RPL_THEMES) # Gamma male popn


#Data set
data <- list(outcome.matrix = outcome.matrix, 
             Pop_m1 = Pop_m1, Pop_m2 = Pop_m2,
             i.spat1 = i.spat1, i.spat2 = i.spat2, 
             i.date1 = i.date1, i.date2 = i.date2
)

## -------------------------------------------------------------------------------------------
### Priors
## From Moraga - gives better estimate of Phi (INLA example gives prec.u 0.2/0.31)
prec_u <- .5/.31
prec_alpha <- .01
phi_u <- .5
phi_alpha <- 2/3

prior_bym2 <- list(
  prec = list(
    prior = "pc.prec",
    param = c(prec_u, prec_alpha)),
  phi = list(
    prior = "pc",
    param = c(phi_u, phi_alpha))
)

prior_bym2 <- list(
  prec = list(
    prior = "pc.prec",
    param = c(prec_u, prec_alpha)),
  phi = list(
    prior = "pc",
    param = c(phi_u, phi_alpha))
)

## From INLA Germany example (time should be a RW)
rw1_u = .2/.31
rw1_alpha = .01
prior_rw1 <- list(theta = list(prior="pc.prec", param=c(rw1_u, rw1_alpha)))

prior_prec <- list(prec = list(prior = "pc.prec",
                               param = c(1, 0.01)))

## -------------------------------------------------------------------------------------------
## Model formula
formula <- outcome.matrix ~ mu -1 +
  Pop_m1 + Pop_m2 + gini1 + gini2 +
  Uninsured1 + Uninsured2 + Production1 + Production2 +
  RUCC1 + RUCC2 + PCP1 + PCP2 + 
  RPL_THEMES1 + RPL_THEMES2 +
  f(i.date1, model = "rw1", hyper = prior_rw1) +
  f(i.date2, copy = "i.date1", fixed = FALSE) +
  f(i.spat1, model = "bym2", graph = g, hyper = prior_bym2) +
  f(i.spat2, copy = "i.spat1", fixed = FALSE)

res1 <- inla(formula, family = c("binomial", "gamma"), data = data, 
             control.compute = list(dic = TRUE), 
             control.predictor = list(compute = TRUE),
             control.fixed = list(expand.factor.strategy = "inla"), 
             control.inla = list(diagonal = 1,
                                 strategy = "adaptive", 
                                 int.strategy = "eb"),
             verbose=TRUE) 

p1 <- plot_fixed_marginals(res1)
ggsave("fixed_hg.pdf", p1)
p1 <- plot_hyper_marginals(res1)
ggsave("hyper_hg.pdf", p1)

time_re1 <- res1$summary.random$i.date1
p1 <- ggplot(time_re1, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = 'gray70') +
  geom_line(aes(y = mean)) + 
  scale_x_continuous("Days since 2020-01-22") + 
  scale_y_continuous("Time random effect") + theme_bw()
print(p1)
time_re2 <- res1$summary.random$i.date2
p2 <- ggplot(time_re2, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = 'gray70') +
  geom_line(aes(y = mean)) + 
  scale_x_continuous("Days since 2020-01-22") + 
  scale_y_continuous("Time random effect") + theme_bw()
print(p2)

ggsave("time_re_hg.pdf", ggarrange(p1, p2))

# dat_sf$u <- res1$summary.random$struct[1:3108, "mean"]

spat_re1 <- res1$summary.random$i.spat1
spat_re2 <- res1$summary.random$i.spat2

out <- list(summary.fixed = res1$summary.fixed, 
            marginals.fixed = res1$marginals.fixed, 
            summary.random = res1$summary.random, 
            marginals.random = res1$marginals.random, 
            summary.hyperpar = res1$summary.hyperpar,
            marginals.hyperpar = res1$marginals.hyperpar)

save(time_re1, time_re2, spat_re1, spat_re2, dat_sf, out,
     file = "./full_model_hg.RData")
## garbage collection
gc()
