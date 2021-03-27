## Preprocessing

## Probably don't need all of these
library(sf)
library(tmap)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(lubridate)
library(spdep)
library(INLA)

## Data first
dat <- read.csv("data/countyTable_timeSeries.csv")
head(dat)

## Shapefile of counties
cnty <- st_read("data/CONUS_counties_1.shp")
# plot(st_geometry(cnty))

# Make up comparable FIPS code for merging
cnty$FIPS = as.numeric(paste0(cnty$STATEFP, cnty$COUNTYFP))
# Index of counties for INLA 
cnty$struct<-1:nrow(cnty)

## Queens adjacency
nb2 <- poly2nb(cnty)
## Convert to INAL graph
nb2INLA("map.adj", nb2)
# g <- inla.read.graph(filename = "map.adj")

## Merge to reorder data to match counties (note that R will sort by default)
dat2 <- merge(cnty, dat, by = "FIPS", sort = FALSE)

## Add date index
dat2$date1 <- as.numeric(ymd(dat2$date) - min(ymd(dat2$date)) + 1)

## Sort by struct and date1
dat2 <- dat2 %>% 
  arrange(struct, date1
          )
## Estimate positivity rate (as percentage)

## First remove zero tests
# dat2 <- dat2 %>%
#   filter(TESTS7AVG_CDC > 0, 
#          CASE7AVG_CDC >= 0)
# dat2$pos_rate <- (dat2$CASE7AVG_CDC / dat2$TESTS7AVG_CDC) * 100
# dat2$pos_rate[dat2$TESTS7AVG_CDC == 0] <- 0

# dat2 <- dat2 %>%
#   filter(!is.na(PERCENTPOS7_CDC))

## Save out
save(dat2, file = "dat2.RData")
