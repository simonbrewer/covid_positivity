## ----setup, include=FALSE-------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## -------------------------------------------------------------------------------------------
library(sf)
library(tmap)
library(dplyr)
library(ggpubr)
library(lubridate)
library(spdep)
library(INLA)
library(INLAutils)


## -------------------------------------------------------------------------------------------
load("dat2.RData")
dim(dat2)

## Convert outcome to proportion
# dat2$TESTSNEW7POS_CDC <- dat2$TESTSNEW7POS_CDC / 100

## -------------------------------------------------------------------------------------------
dat.sub <- dat2 %>%
  filter(date == "2021-01-01") 


## -------------------------------------------------------------------------------------------
tm_shape(dat.sub) +
  tm_fill("TESTSNEW7POS_CDC")


## -------------------------------------------------------------------------------------------
g <- inla.read.graph(filename = "./map.adj")


## -------------------------------------------------------------------------------------------
### Priors
## From Moraga - gives better estimate of Phi
u <- .5/.31
alpha <- .01
phi.u <- .5
phi.alpha <- 2/3

prior_bym2 <- list(
  prec = list(
    prior = "pc.prec",
    param = c(u, alpha)),
  phi = list(
    prior = "pc",
    param = c(phi.u, phi.alpha))
)

## -------------------------------------------------------------------------------------------
f0 <- TESTSNEW7POS_CDC ~
  1 +
  f(state, model = "iid") 

f0 <- test2 ~
  1 +
  f(state, model = "iid") 

## -------------------------------------------------------------------------------------------
res0 <- inla(f0, data = dat.sub, family = 'binomial',
             control.compute = list(dic = TRUE, waic = TRUE), 
             #Ntrials = ceiling(dat.sub$TESTSNEW7AVG),
             verbose =TRUE
)


## -------------------------------------------------------------------------------------------
summary(res0)

## -------------------------------------------------------------------------------------------
f1 <- TESTSNEW7POS_CDC ~
  1 +
  f(struct, model = "bym2", graph = g, hyper = prior_bym2) 

# dat.sub <- dat.sub %>% 
#   filter(TESTSNEW7AVG > 0)

## -------------------------------------------------------------------------------------------
res1 <- inla(f1, data = dat.sub, family = 'binomial',
             control.family=list(link='logit'),
             control.predictor = list(compute = TRUE),
             control.compute = list(dic = TRUE, waic = TRUE), 
             Ntrials = dat.sub$TESTSNEW7AVG,
             verbose =TRUE
)


## -------------------------------------------------------------------------------------------
summary(res1)


## -------------------------------------------------------------------------------------------
p <- autoplot(res1)
print(p)

## -------------------------------------------------------------------------------------------
res1$summary.random$struct


## -------------------------------------------------------------------------------------------
dat.sub$v <- res1$summary.random$struct[1:3108, "mean"]
dat.sub$u <- res1$summary.random$struct[3109:6216, "mean"]


## -------------------------------------------------------------------------------------------
tm_shape(dat.sub) +
  tm_fill("u")


## -------------------------------------------------------------------------------------------
f2 <- TESTSNEW7POS_CDC ~
  Pop_m + gini + Uninsured + Production + RUCC + PCP + nursing + universities + RPL_THEMES + 
  f(state, model = "iid", 
    constr=TRUE) 


## -------------------------------------------------------------------------------------------
res2 <- inla(f2, data = dat.sub,
             control.predictor = list(compute = TRUE),
             control.compute = list(dic = TRUE, waic = TRUE), verbose = TRUE, 
             control.inla = list(strategy = "simplified.laplace", 
                                 int.strategy="eb")
)

## -------------------------------------------------------------------------------------------
summary(res2)
plot_fixed_marginals(res2)
plot_hyper_marginals(res2)
plot_random_effects(res2)
plot_random_effects(res2, type = 'boxplot')

## -------------------------------------------------------------------------------------------
## BYM2 model
f3 <- TESTSNEW7POS_CDC ~
  Pop_m + gini + Uninsured + Production + RUCC + PCP + nursing + universities + RPL_THEMES + 
  f(struct, model = "bym2", graph = g, hyper = prior_bym2, 
    constr=TRUE) 

## -------------------------------------------------------------------------------------------
res3 <- inla(f3, data = dat.sub,
             control.predictor = list(compute = TRUE),
             control.compute = list(dic = TRUE, waic = TRUE), verbose = TRUE, 
             control.family = list( 
               hyper = list( 
                 prec = list( 
                   initial = -1))),
             control.inla = list(strategy = "simplified.laplace", 
                                 int.strategy="eb")
)

## -------------------------------------------------------------------------------------------
summary(res3)
plot_fixed_marginals(res3)
plot_hyper_marginals(res3)


## -------------------------------------------------------------------------------------------
## BYM2 model + state effect
f4 <- TESTSNEW7POS_CDC ~
  Pop_m + gini + Uninsured + Production + RUCC + PCP + nursing + universities + RPL_THEMES + 
  f(state, model = "iid") +
  f(struct, model = "bym2", graph = g, hyper = prior_bym2) 

## -------------------------------------------------------------------------------------------
res4 <- inla(f4, data = dat.sub,
             control.predictor = list(compute = TRUE),
             control.compute = list(dic = TRUE, waic = TRUE), verbose = TRUE, 
             control.family = list( 
               hyper = list( 
                 prec = list( 
                   initial = -1))),
             control.inla = list(strategy = "simplified.laplace", 
                                 int.strategy="eb")
)


## -------------------------------------------------------------------------------------------
summary(res4)
plot_fixed_marginals(res4)
plot_hyper_marginals(res4)
plot_random_effects(res4, type = 'boxplot')

stop()
## Testing below here

### Priors
# prior_bym2 <- list(
#   prec = list(
#     prior = "pc.prec",
#     param = c(0.5 / 0.31, 0.01)),
#   phi = list(
#     prior = "pc",
#     param = c(0.5, 2 / 3))
# )

## From INLA help
prec.u <- .2/.31
prec.alpha <- .01
phi.u <- .5
phi.alpha <- 2/3

## From Moraga - gives better estimate of Phi
# u <- .5/.31
# alpha <- .01
# phi.u <- .5
# phi.alpha <- 2/3

f3.1 <- PERCENTPOS7_CDC ~
  Pop_o_60 + Pop_m + Pop_black + Pop_white + Pop_AmIndAlNat + Pop_asia + Pop_NaHaPaIs +
  Income + Bachelor + Disabled + Unemployed + Uninsured + Poverty + Production +
  RUCC + hospitals + nursing + universities +
  f(state, model = "iid") +
  f(struct, model = "bym2", graph = g, 
    constr=TRUE, scale.model=TRUE, 
    hyper=list(phi=list(prior='pc',
                        param=c(phi.u, phi.alpha),
                        initial=-3),
               prec=list(prior='pc.prec',
                         param=c(prec.u, prec.alpha),
                         initial=5))) 

## -------------------------------------------------------------------------------------------
res3.1 <- inla(f3.1, data = dat.sub,
               control.predictor = list(compute = TRUE),
               control.compute = list(dic = TRUE, waic = TRUE), verbose = TRUE, 
               control.family = list( 
                 hyper = list( 
                   prec = list( 
                     initial = -1))),
               control.inla = list(strategy = "simplified.laplace", 
                                   int.strategy="eb")
)
summary(res3.1)
plot_fixed_marginals(res3.1)
plot_hyper_marginals(res3.1)
# formula2 <- agr ~ age + sex + age:sex +
#   f(geo, model='bym2', graph=border,
#     constr=TRUE, scale.model=TRUE,
#     hyper=list(phi=list(prior='pc',
#                         param=c(phi.u, phi.alpha),
#                         initial=-3),
#                prec=list(prior='pc.prec',
#                          param=c(u, alpha),
#                          initial=5)))
