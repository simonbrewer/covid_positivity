## ----setup, include=FALSE-------------------------------------------------------------------

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

prior_prec <- list(prec = list(prior = "pc.prec",
                               param = c(1, 0.01)))

## -------------------------------------------------------------------------------------------
## Model 1: just the intercept
f1 <- TESTSNEW7POS_CDC ~ scale(Pop_m) + scale(gini) + scale(Uninsured) + 
  scale(Production) + RUCC + 
  scale(PCP) + scale(RPL_THEMES) +
  f(date1, model = "rw1", hyper = prior_rw1, 
    scale.model = TRUE) +
  f(struct, model = "bym2", graph = g, hyper = prior_bym2, 
    scale.model = TRUE)

dat2$TESTSNEW7POS_CDC <- dat2$TESTSNEW7POS_CDC / 100
cens <- min(dat2$TESTSNEW7POS_CDC[dat2$TESTSNEW7POS_CDC > 0])

## diagonal = 0
system.time(
  res1 <- inla(f1, data = dat2,
               family = "beta",
               control.family=list(link='logit', 
                                   beta.censor.value = cens),
               #Ntrials = dat2$TESTSNEW7AVG,
               control.predictor = list(compute = TRUE),
               control.compute = list(dic = TRUE, waic = TRUE, 
                                      openmp.strategy = "pardiso"),
               control.fixed = list(prec.intercept = 1),
               #control.mode = list(result = res1, restart = TRUE),
               control.inla = list(int.strategy = "eb",
                                   diagonal = 0),
               verbose = TRUE, num.threads = 4:-1,
  )
)

p1 <- plot_fixed_marginals(res1)
ggsave("fixed_cb.pdf", p1)
p1 <- plot_hyper_marginals(res1)
ggsave("hyper_cb.pdf", p1)

time.re <- res1$summary.random$date1
p1 <- ggplot(time.re, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = 'gray70') +
  geom_line(aes(y = mean)) + 
  scale_x_continuous("Days since 2020-01-22") + 
  scale_y_continuous("Time random effect") + theme_bw()
print(p1)
ggsave("time_re_cb.pdf", p1)

## Crude trend plot in odds
time.re$mean_odds <- exp(time.re$mean)
time.re$cilo_odds <- exp(time.re$`0.025quant`)
time.re$cihi_odds <- exp(time.re$`0.975quant`)
p2 <- ggplot(time.re, aes(x = ID)) + 
  geom_ribbon(aes(ymin = cilo_odds, ymax = cihi_odds), fill = 'gray70') +
  geom_line(aes(y = mean_odds)) + 
  scale_x_continuous("Days since 2020-01-22") + 
  scale_y_continuous("Time random effect") + theme_bw()
ggsave("time_re_odds_cb.pdf", p2)


dat_sf$u <- res1$summary.random$struct[1:3108, "mean"]

out <- list(summary.fixed = res1$summary.fixed, 
            marginals.fixed = res1$marginals.fixed, 
            summary.random = res1$summary.random, 
            marginals.random = res1$marginals.random, 
            summary.hyperpar = res1$summary.hyperpar,
            marginals.hyperpar = res1$marginals.hyperpar)

save(time.re, dat_sf, out,
     file = "./full_model_cb.RData")
## garbage collection
gc()
## -------------------------------------------------------------------------------------------

stop()
## -------------------------------------------------------------------------------------------
## Model 2: include state effect
f2 <- TESTSNEW7POS_CDC ~ scale(Pop_m) + scale(gini) + scale(Uninsured) + 
  scale(Production) + RUCC + 
  scale(PCP) + scale(RPL_THEMES) +
  f(date1, model = "rw1", hyper = prior_rw1, 
    scale.model = TRUE) +
  f(struct, model = "bym2", graph = g, hyper = prior_bym2, 
    scale.model = TRUE) + 
  f(state, model = "iid", constr = TRUE)

dat2$TESTSNEW7POS_CDC <- dat2$TESTSNEW7POS_CDC / 100
cens <- min(dat2$TESTSNEW7POS_CDC[dat2$TESTSNEW7POS_CDC > 0])

## diagonal = 0
system.time(
  res2 <- inla(f2, data = dat2,
               family = "beta",
               control.family=list(link='logit', 
                                   beta.censor.value = cens),
               #Ntrials = dat2$TESTSNEW7AVG,
               control.predictor = list(compute = TRUE),
               control.compute = list(dic = TRUE, waic = TRUE, 
                                      openmp.strategy = "pardiso"),
               control.fixed = list(prec.intercept = 1),
               #control.mode = list(result = res1, restart = TRUE),
               control.inla = list(int.strategy = "eb",
                                   diagonal = 0),
               verbose = TRUE, num.threads = 4:-1,
  )
)

p1 <- plot_fixed_marginals(res2)
ggsave("fixed_cb_s.pdf", p1)
p1 <- plot_hyper_marginals(res2)
ggsave("hyper_cb_s.pdf", p1)

time.re <- res2$summary.random$date1
p1 <- ggplot(time.re, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = 'gray70') +
  geom_line(aes(y = mean)) + 
  scale_x_continuous("Days since 2020-01-22") + 
  scale_y_continuous("Time random effect") + theme_bw()
print(p1)
ggsave("time_re_cb_s.pdf", p1)

## Crude trend plot in odds
time.re$mean_odds <- exp(time.re$mean)
time.re$cilo_odds <- exp(time.re$`0.025quant`)
time.re$cihi_odds <- exp(time.re$`0.975quant`)
p2 <- ggplot(time.re, aes(x = ID)) + 
  geom_ribbon(aes(ymin = cilo_odds, ymax = cihi_odds), fill = 'gray70') +
  geom_line(aes(y = mean_odds)) + 
  scale_x_continuous("Days since 2020-01-22") + 
  scale_y_continuous("Time random effect") + theme_bw()
ggsave("time_re_odds_cb_s.pdf", p2)


dat_sf$u <- res2$summary.random$struct[1:3108, "mean"]

out <- list(summary.fixed = res2$summary.fixed, 
            marginals.fixed = res2$marginals.fixed, 
            summary.random = res2$summary.random, 
            marginals.random = res2$marginals.random, 
            summary.hyperpar = res2$summary.hyperpar,
            marginals.hyperpar = res2$marginals.hyperpar)

save(time.re, dat_sf, out,
     file = "./full_model_cb_s.RData")
## garbage collection
gc()
## -------------------------------------------------------------------------------------------

