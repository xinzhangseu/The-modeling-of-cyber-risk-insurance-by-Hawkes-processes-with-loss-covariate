##This code is used to give the detailed description of the data sets. 
rm(list=ls())
library(ggplot2)
#library(tidyverse)
library(lubridate)
library(dplyr)
setwd("E:/paper/rdata/nonstationary")
load(file="E:/paper/rdata/nonstationary/prc_hack_15_bus.rda")
load(file="E:/paper/rdata/nonstationary/prc_hack_15_med.rda")
save(prc_hack_bus_arr_interval,file="E:/paper/rdata/nonstationary/prc_hack_bus_arr_interval.rda")
save(prc_hack_med_arr,file="E:/paper/rdata/nonstationary/prc_hack_med_arr.rda")
prc_hack_2015=filter(prc_hack,Year.of.Breach=="2015"|
                       Year.of.Breach=="2016"|
                       Year.of.Breach=="2017"|
                       Year.of.Breach=="2018")

prc_hack_15_bus=filter(prc_hack_2015,Type.of.organization=="BSF"|
                      Type.of.organization=="BSO"|
                      Type.of.organization=="BSF")
prc_hack_15_med=filter(prc_hack_2015,Type.of.organization=="MED")
###### statistical description
plot(x = prc_hack_15_med$Date.Made.Public, y = prc_hack_15_med$Total.Records,type="h",
     xlab="Time of release",ylab="Number of records",
     cex.lab=1.5,cex.axis=1.5, cex.main=1.5)  
plot(x = prc_hack_15_bus$Date.Made.Public, y = prc_hack_15_bus$Total.Records,type="h",
     xlab="Time of release",ylab="Number of records",
     cex.lab=1.5,cex.axis=1.5, cex.main=1.5) 
#####interval time
prc_hack_bus_arr=arrange(prc_hack_15_bus,prc_hack_15_bus[,1])
prc_hack_med_arr=arrange(prc_hack_15_med,prc_hack_15_med[,1])
prc_hack_med_arr_interval=diff(prc_hack_med_arr$Date.Made.Public,lag=1)
prc_hack_bus_arr_interval=diff(prc_hack_bus_arr$Date.Made.Public,lag=1)
######QQ plot
p <- ppoints(525)
q <- quantile(prc_hack_med_arr_interval,p=p)
plot(qexp(p),q, type="p",main="Q-Q Plot",
     xlab="Theoretical Quantiles",ylab="Sample Quantiles",
     cex.lab=1.5,cex.axis=1.5, cex.main=1.5)  
qqline(q, distribution=qexp,col="red", lty=2)

p <- ppoints(126)
q <- quantile(prc_hack_bus_arr_interval,p=p)
plot(qexp(p),q, type="p",main="Q-Q Plot",
     xlab="Theoretical Quantiles",ylab="Sample Quantiles",
     cex.lab=1.5,cex.axis=1.5, cex.main=1.5)  
qqline(q, distribution=qexp,col="red", lty=2)
######
acf(prc_hack_med_arr_interval,lag=10,pl=TRUE,main="",
    cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
acf(prc_hack_bus_arr_interval,lag=10,pl=TRUE,main="",
    cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
#####
hist(as.numeric(prc_hack_med_arr_interval),xlab="Inter-arrival times",
     ylab="Frequency",main="",
     cex.lab=1.5,cex.axis=1.5, cex.main=1.5)

hist(as.numeric(prc_hack_bus_arr_interval),xlab="Inter-arrival times",
     ylab="Frequency",main="",
     cex.lab=1.5,cex.axis=1.5, cex.main=1.5)