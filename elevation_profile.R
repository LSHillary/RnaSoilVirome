rm(list=ls())
setwd("~/Dropbox/PhD/Chapter 4 - Soil RNA Virome/R")

library(tidyverse)  # ggplot() fortify()
library(data.table)

elevation_data <- fread("elevation profile.csv", header = T)
ggplot(elevation_data, aes(x = Position, y = Elevation)) +
  geom_area(fill = "chartreuse 4") +
  geom_point(data = elevation_data[elevation_data$V3 != "",], aes(x = Position, y = Elevation, shape = V3), colour = "black", size = 3) +
  geom_text(data = elevation_data[elevation_data$V3 != "",], aes(x = Position, y = Elevation, label = V3), colour = "black", size = 3, vjust = -1, hjust = 0) +
  scale_shape_manual(values=c(16, 16, 16, 16, 16))+
  ylim(0, 500) +
  xlim(0,3500) +
  theme_bw() +
  theme(legend.position = "none") +
  labs (x = "Distance (m)",
        y = "Elevation (m)")
