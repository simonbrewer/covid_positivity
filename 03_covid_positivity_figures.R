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
load("./full_model.RData")

p1 <- ggplot(time.re, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = 'gray70') +
  geom_line(aes(y = mean)) + 
  scale_x_continuous("Days since 2020-01-22") + 
  scale_y_continuous("Smoother") + theme_bw() +
  ggtitle("Day random effect")
ggsave(p1, file = "time_random_effects.pdf")
print(p1)

p2 <- ggplot(time.re, aes(x = ID)) + 
  geom_ribbon(aes(ymin = cilo_odds, ymax = cihi_odds), fill = 'gray70') +
  geom_line(aes(y = mean_odds)) + 
  scale_x_continuous("Days since 2020-01-22") + 
  scale_y_continuous("Smoother") + theme_bw() +
  ggtitle("Day random effect")
# ggsave(p1, file = "time_random_effects.pdf")
print(p2)

m1 <- tm_shape(dat_sf) + tm_fill("u", n = 9) +
  tm_layout(main.title = "Spatial RE")
print(m1)

tmap_save(m1, file = "spatial_random_effects.pdf")

## -------------------------------------------------------------------------------------------
