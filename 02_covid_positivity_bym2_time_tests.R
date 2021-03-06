## ----setup, include=FALSE-------------------------------------------------------------------

## -------------------------------------------------------------------------------------------
#library(sf)
#library(tmap)
library(dplyr)
library(ggpubr)
library(lubridate)
library(spdep)
library(INLA)
library(INLAutils)
library(ggregplot)

# inla.setOption(pardiso.license = "~/pardiso.lic")

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

## From INLA Germany example (time should be a RW)
rw1_u = .2/.31
rw1_alpha = .01
prior_rw1 <- list(theta = list(prior="pc.prec", param=c(rw1_u, rw1_alpha)))

## -------------------------------------------------------------------------------------------
## Subset data
# dat.sub <- dat2 %>%
#   filter(date1 < 150) %>%
#   # select(struct, date1, PERCENTPOS7_CDC) %>%
#   mutate(lppos7 = log(PERCENTPOS7_CDC+1e-1))
dat.sub <- dat2 %>%
  filter(date1 < 100)
# dat.sub <- dat2 %>%
#   filter(date1 > 100 & date1 < 150) %>%
#   select(TESTSNEW7POS_CDC, Pop_m, gini, Uninsured, Production,
#          RUCC, PCP, nursing, universities, RPL_THEMES, struct, date1)

prior_prec <- list(prec = list(prior = "pc.prec",
                               param = c(1, 0.01)))

## -------------------------------------------------------------------------------------------
## Model 1: just the intercept
f1 <- TESTSNEW7POS_CDC ~ 1

## Model 2: covariates
f2 <- TESTSNEW7POS_CDC ~ Pop_m + gini + Uninsured + Production + RUCC + 
  PCP + nursing + universities + RPL_THEMES

## Model 3: intercept + time
f3 <- TESTSNEW7POS_CDC ~ 1 +
  f(date1, model = "rw1", hyper = prior_rw1)

## Model 4: covariates + time
f4 <- TESTSNEW7POS_CDC ~ Pop_m + gini + Uninsured + Production + RUCC + 
  PCP + nursing + universities + RPL_THEMES +
  f(date1, model = "rw1", hyper = prior_rw1)

## Model 5: intercept + space
f5 <- TESTSNEW7POS_CDC ~ 1 +
  f(struct, model = "bym2", graph = g, hyper = prior_bym2, 
    scale.model = TRUE, diagonal = 1e-2)

## Model 6: covariates + space
f6 <- TESTSNEW7POS_CDC ~ Pop_m + gini + Uninsured + Production + RUCC + 
  PCP + nursing + universities + RPL_THEMES +
  f(struct, model = "bym2", graph = g, hyper = prior_bym2, 
    scale.model = TRUE, diagonal = 1e-1)

## Model 7: intercept + space and time
f7 <- TESTSNEW7POS_CDC ~ 1 +
  f(date1, model = "rw1", hyper = prior_rw1, 
    scale.model = TRUE, diagonal = 1e-1) +
  f(struct, model = "bym2", graph = g, hyper = prior_bym2, 
    scale.model = TRUE, diagonal = 1e-1)

# f8 <- TESTSNEW7POS_CDC ~ Pop_m + gini + Uninsured + Production + RUCC + 
#   PCP + nursing + universities + RPL_THEMES
## Model 8: intercept + space and time
f8 <- TESTSNEW7POS_CDC ~ Pop_m + gini + Uninsured + Production + RUCC + 
  PCP + RPL_THEMES +
  f(date1, model = "rw1", hyper = prior_rw1, 
    scale.model = TRUE, diagonal = 5e-1) +
  f(struct, model = "bym2", graph = g, hyper = prior_bym2, 
    scale.model = TRUE, diagonal = 5e-1)

## Model 8: intercept + space and time (scaled)
## Added prior on intercept
f8 <- TESTSNEW7POS_CDC ~ scale(Pop_m) + scale(gini) + scale(Uninsured) + 
  scale(Production) + RUCC + 
  scale(PCP) + scale(RPL_THEMES) +
  f(date1, model = "rw1", hyper = prior_rw1, 
    scale.model = TRUE, diagonal = 1e-2) +
  f(struct, model = "bym2", graph = g, hyper = prior_bym2, 
    scale.model = TRUE, diagonal = 1e-2)

## f9 - constraint to zero
## constr = TRUE
system.time(
  res1 <- inla(f8, data = dat2,
               control.predictor = list(compute = TRUE),
               control.compute = list(dic = TRUE, waic = TRUE), 
               control.fixed = list(prec.intercept = 0.1),
               control.inla = list(strategy = "adaptive",
                                   int.strategy="eb"),
               verbose = TRUE, num.threads = 32,
  )
)

plot_fixed_marginals(res1)
plot_hyper_marginals(res1)

time.re <- res1$summary.random$date1
p1 <- ggplot(time.re, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = 'gray70') +
  geom_line(aes(y = mean)) + 
  scale_x_continuous("Days since 2020-06-02") + 
  scale_y_continuous("Time random effect") + theme_bw()
print(p1)


## garbage collection
gc()
## -------------------------------------------------------------------------------------------

