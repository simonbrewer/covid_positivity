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
## County IID
# dat.sub <- dat2 %>%
#   filter(date1 < 150) %>%
#   # select(struct, date1, PERCENTPOS7_CDC) %>%
#   mutate(lppos7 = log(PERCENTPOS7_CDC+1e-1))
dat.sub <- dat2 %>%
  filter(date1 < 200)
dat.sub <- dat2 %>%
  filter(date1 > 100 & date1 < 250) %>%
  select(TESTSNEW7POS_CDC, Pop_m, gini, Uninsured, Production,
         RUCC, PCP, nursing, universities, RPL_THEMES, struct, date1)

prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(1, 0.01)))
# prior.prec <- list(prec = list(initial = log(0.000001),
#                                fixed = TRUE))

f1 <- TESTSNEW7POS_CDC ~ 
  1 +
  f(struct, model = "bym2", graph = g, hyper = prior_bym2) +
  f(date1, model = "rw1", hyper = prior_rw1)

# system.time(
#   res1 <- inla(f1, data = dat.sub,
#                control.predictor = list(compute = TRUE),
#                control.compute = list(dic = TRUE, waic = TRUE),
#                # control.fixed = list(prec.intercept = 0.1),
#                control.inla = list(diagonal = 0,
#                                    strategy = "simplified.laplace", 
#                                    int.strategy="eb"), 
#                verbose = TRUE, num.threads = 4,
#   )
# )
# 
# summary(res1)
# # plot_fixed_marginals(res1)
# plot_hyper_marginals(res1)
# Efxplot(res1)
# 
# time.re <- res1$summary.random$date1
# p1 <- ggplot(time.re, aes(x = ID)) + 
#   geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = 'gray70') +
#   geom_line(aes(y = mean)) + 
#   scale_x_continuous("Days since 2020-06-02") + 
#   scale_y_continuous("Time random effect") + theme_bw()
# print(p1)
# 
# dat_sf$u <- res1$summary.random$struct[1:3108, "mean"]
# tm_shape(dat_sf) + tm_fill("u")
# save(res1, file = "output/null_model.RData")
# rm(res1)
## garbage collection
gc()
## -------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------
## Full model
f2 <- TESTSNEW7POS_CDC ~
  Pop_m + gini + Uninsured + Production + RUCC + 
  PCP + nursing + universities + RPL_THEMES + 
  f(struct, model = "bym2", graph = g, hyper = prior_bym2, 
    scale.model = TRUE, diagonal = 1e-2) +
  f(date1, model = "rw1", hyper = prior_rw1, 
    scale.model = TRUE, diagonal = 1e-2)

system.time(
  res2 <- inla(f2, data = dat2,
               control.predictor = list(compute = TRUE),
               control.compute = list(dic = TRUE, waic = TRUE,
                                      openmp.strategy = "pardiso.serial"),
               # control.fixed = list(prec.intercept = 0.1),
               control.inla = list(diagonal = 0,
                                   tolerance = 1e-6,
                                   strategy = "simplified.laplace", 
                                   int.strategy="eb"), 
               verbose = TRUE, num.threads = 2
  )
)

summary(res2)
# plot_fixed_marginals(res1)
pdf("fixed_effects.pdf")
plot_hyper_marginals(res2)
Efxplot(res2)
dev.off()

time.re <- res2$summary.random$date1
p1 <- ggplot(time.re, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = 'gray70') +
  geom_line(aes(y = mean)) + 
  scale_x_continuous("Days since 2020-06-02") + 
  scale_y_continuous("Time random effect") + theme_bw()
print(p1)

dat_sf$u <- res2$summary.random$struct[1:3108, "mean"]
tm_shape(dat_sf) + tm_fill("u", palette = "-magma")

save(time.re, dat_sf, file = "output/full_model.RData")
# save(res2, file = "output/full_model.RData")
## garbage collection
gc()
## -------------------------------------------------------------------------------------------
