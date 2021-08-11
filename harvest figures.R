

total_harvest_t <- c(49467, 44362, 66970, 63430, 50352, 94849)
total_harvest_c <- c(16401, 7829, 47140, 40373, 17515, 21717)


total_harvest <- total_harvest_t + total_harvest_c
total_harvest


th <- c(924392, 603265, 816261, NA, 625274, 1233189, NA, NA, 664101, 501424, 669282, NA, NA, 619356, NA, NA, NA, NA, 286856, 172909, 222383, NA, NA, NA, 186927, 158372, NA, NA, NA, 85420, 114311, 119427, 66585, 53283)

th <- c(th, total_harvest)

length(th) + 1976

year <- 1977:2016

length(th)
length(year)

plot(x = year, y = th/1000, 
     ylab = "Total murre harvest (x1000)", xlab = "Year",
     type = "h", lwd = 5)
abline(h=seq(100, 1200, 100), lty = 3)
abline(h=375, lty = 1, col = "red", lwd = 17)
abline(v=1993, lty =2)


##############################################
#COMU trend figure
###############################################

library(tidyr)
library(dplyr)
library(plyr)

setwd("C:\\Users\\robertsong\\Documents\\seabirds\\murres\\harvest\\2015 assessment")

COMU <- read.csv("COMU trends NL.csv", header = TRUE)

COMU

COMU <- COMU %>% gather(COMU, "Count", -Year)

names(COMU)[2] <- "Colony"

COMU
library(lattice)

xyplot(log10(Count) ~ Year|Colony, data = subset(COMU, Year > 1971 & Colony != "The.Doughboy"), type = c("p"),
       ylab = "Number of pairs",
       scales=list(y=list(labels=c(1, 10, 100, 1000, 10000, "100000"))))

str(COMU)


