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
library(ggregplot)

## -------------------------------------------------------------------------------------------
# load("dat2.RData")
# ## Time index
# dat2$date1 <- as.numeric(ymd(dat2$date) - min(ymd(dat2$date)) + 1)
## SF object for plotting
dat_sf <- dat2 %>%
  filter(date1 == 1)

load("./output/full_model.RData")

p1 <- ggplot(time.re, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = 'gray70') +
  geom_line(aes(y = mean)) + 
  scale_x_continuous("Days since 2020-06-02") + 
  scale_y_continuous("Smoother") + theme_bw() +
  ggtitle("Day random effect")
ggsave(p1, file = "time_random_effects.pdf")
print(p1)

m1 <- tm_shape(dat_sf) + tm_fill("u", palette = "-magma") +
  tm_layout(main.title = "Spatial RE")
print(m1)

dat_sf$u_sig <- rep(0, nrow(dat_sf))
dat_sf$u_sig[dat_sf$u_lo > 0] <- 1
dat_sf$u_sig[dat_sf$u_hi < 0] <- -1
dat_sf$u_sig <- as.factor(dat_sf$u_sig)
m2 <- tm_shape(dat_sf) + tm_fill("u_sig", palette = "Set3") +
  tm_layout(main.title = "Spatial RE")
print(m2)

tmap_save(tmap_arrange(m1, m2), file = "spatial_random_effects.pdf")

## -------------------------------------------------------------------------------------------
