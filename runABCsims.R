###############################################################
## Program to run the country simulations using ABC approach ##
###############################################################

source("utils.R")
library(readxl)
library(ggplot2)
library(reshape)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)
library(grid)
library(scales)
library(dplyr)
library(data.table)
library(xlsx)
library(kableExtra)
library(tidyr)

# Read in the data #
####################
# prevalence survey estimates, with years and CIs
dat.prev <- read_excel("Data/PrevalenceSurveyData.xlsx",sheet="Sheet1")

# data on treatments 
dat.cure <- read_excel("Data/AnnualTreatments.xlsx",sheet="Sheet1")

# number of simulations to run
M <- 1000


#############
## VIETNAM ##
#############
country.tmp <- "Vietnam"
viet.ann.rate <- 1-0.045
tmp.pop <- read.csv("Data/vietnam-pop.csv")
y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]

col.names <- c("Duration of treated disease","LCL","UCL","Assumed duration of untreated disease")

# run model for untreated TB assumed to be between 3 and 7 years
sim.dat.viet <- list()
dur.ests.viet<- data.frame(matrix(0,nr=length(c(3:7)),nc=4))
prev.ests.1 <- prev.ests.2 <- data.frame(matrix(0,nr=M*length(c(3:7)),nc=y2-y1+3))
for(j in 3:7){
  tmp.dat <- getDurationABC(country=country.tmp,pop.data=tmp.pop,duration.range=seq(1,6,0.01),unDxTBdur=j,M=M,ABC.cut.off=0.15)
  sim.dat.viet[[j-2]] <- tmp.dat
  
  # add estimates of duration
  dur.ests.viet[j-2,] <- c(tmp.dat$duration,tmp.dat$CIs,tmp.dat$Duration.u)
  
  # get prevalence estimates from each simulation
  tmp.prev <- getSummaryData(getDurationABC.output=tmp.dat)
  
  prev.ests.1[((j-3)*M+1):((j-3)*M+M),] <- cbind(tmp.prev$ann.prev.1,rep(j,M),rep(tmp.dat$duration,M))
  prev.ests.2[((j-3)*M+1):((j-3)*M+M),] <- cbind(tmp.prev$ann.prev.2,rep(j,M),rep(tmp.dat$duration,M))
  
}

names(dur.ests.viet) <- col.names
names(prev.ests.1) <- names(prev.ests.2) <- c(as.character(seq(y1,y2,1)),"Assumed duration of untreated disease","Duration of treated disease")
prev.ests.1.v <- prev.ests.1
prev.ests.2.v <- prev.ests.2

save.image()

###############
## INDONESIA ##
###############
country.tmp <- "Indonesia"
indo.ann.rate <- getRateofChange(prev1=120,prev2=180,year1=2004,year2=2013.5)
tmp.pop <- read.csv("Data/indo-pop.csv")

y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]

# run model for untreated TB assumed to be between 3 and 7 years
sim.dat.indo <- list()
dur.ests.indo<- data.frame(matrix(0,nr=length(c(3:7)),nc=4))
prev.ests.1 <- prev.ests.2 <- data.frame(matrix(0,nr=M*length(c(3:7)),nc=y2-y1+3))
for(j in 3:7){
  tmp.dat <- getDurationABC(country=country.tmp,pop.data=tmp.pop,duration.range=seq(1,6,0.01),unDxTBdur=j,M=M,ABC.cut.off=0.15)
  sim.dat.indo[[j-2]] <- tmp.dat
  
  # add estimates of duration
  dur.ests.indo[j-2,] <- c(tmp.dat$duration,tmp.dat$CIs,tmp.dat$Duration.u)
  
  # get prevalence estimates from each simulation
  tmp.prev <- getSummaryData(getDurationABC.output=tmp.dat)
  
  prev.ests.1[((j-3)*M+1):((j-3)*M+M),] <- cbind(tmp.prev$ann.prev.1,rep(j,M),rep(tmp.dat$duration,M))
  prev.ests.2[((j-3)*M+1):((j-3)*M+M),] <- cbind(tmp.prev$ann.prev.2,rep(j,M),rep(tmp.dat$duration,M))
  
}

names(dur.ests.indo) <- col.names
names(prev.ests.1) <- names(prev.ests.2) <- c(as.character(seq(y1,y2,1)),"Assumed duration of untreated disease","Duration of treated disease")
prev.ests.1.i <- prev.ests.1
prev.ests.2.i <- prev.ests.2

#################
## PHILIPPINES ##
#################
country.tmp <- "Philippines"
rate.change.10y <- 0.105
phil.prev1 <- dat.prev$`1st Prevalence`[dat.prev$Country==country.tmp] 
phil.ann.rate <- getRateofChange(prev1=phil.prev1,prev2=982,year1=2007,year2=2016)

tmp.pop <- read.csv("Data/phil-pop.csv")

y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]

# run model for untreated TB assumed to be between 3 and 7 years
sim.dat.phil <- list()
dur.ests.phil <- data.frame(matrix(0,nr=length(c(3:7)),nc=4))
prev.ests.1 <- prev.ests.2 <- data.frame(matrix(0,nr=M*length(c(3:7)),nc=y2-y1+3))
for(j in 3:7){
  tmp.dat <- getDurationABC(country=country.tmp,pop.data=tmp.pop,duration.range=seq(1,6,0.01),unDxTBdur=j,M=M,ABC.cut.off=0.15)
  sim.dat.phil[[j-2]] <- tmp.dat
  
  # add estimates of duration
  dur.ests.phil[j-2,] <- c(tmp.dat$duration,tmp.dat$CIs,tmp.dat$Duration.u)
  
  # get prevalence estimates from each simulation
  tmp.prev <- getSummaryData(getDurationABC.output=tmp.dat)
  
  prev.ests.1[((j-3)*M+1):((j-3)*M+M),] <- cbind(tmp.prev$ann.prev.1,rep(j,M),rep(tmp.dat$duration,M))
  prev.ests.2[((j-3)*M+1):((j-3)*M+M),] <- cbind(tmp.prev$ann.prev.2,rep(j,M),rep(tmp.dat$duration,M))
  
}

names(dur.ests.phil) <- col.names
names(prev.ests.1) <- names(prev.ests.2) <- c(as.character(seq(y1,y2,1)),"Assumed duration of untreated disease","Duration of treated disease")
prev.ests.1.p <- prev.ests.1
prev.ests.2.p <- prev.ests.2

save.image()

#################
##    CHINA    ##
#################
country.tmp <- "China"
china.perc.change <- 0.35
china.prev1 <- dat.prev$`1st Prevalence`[dat.prev$Country==country.tmp]
china.ann.rate <- getRateofChange(prev1=china.prev1,prev2=108,year1=2000,year2=2010)

tmp.pop <- read.csv("Data/china-pop.csv")


y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]

# run model for untreated TB assumed to be between 3 and 7 years
sim.dat.china <- list()
dur.ests.china <- data.frame(matrix(0,nr=length(c(3:7)),nc=4))
prev.ests.1 <- prev.ests.2 <- data.frame(matrix(0,nr=M*length(c(3:7)),nc=y2-y1+3))
for(j in 3:7){
  tmp.dat <- getDurationABC(country=country.tmp,pop.data=tmp.pop,duration.range=seq(1,6,0.01),unDxTBdur=j,M=M,ABC.cut.off=0.15)
  sim.dat.china[[j-2]] <- tmp.dat
  
  # add estimates of duration
  dur.ests.china[j-2,] <- c(tmp.dat$duration,tmp.dat$CIs,tmp.dat$Duration.u)
  
  # get prevalence estimates from each simulation
  tmp.prev <- getSummaryData(getDurationABC.output=tmp.dat)
  
  prev.ests.1[((j-3)*M+1):((j-3)*M+M),] <- cbind(tmp.prev$ann.prev.1,rep(j,M),rep(tmp.dat$duration,M))
  prev.ests.2[((j-3)*M+1):((j-3)*M+M),] <- cbind(tmp.prev$ann.prev.2,rep(j,M),rep(tmp.dat$duration,M))
  
}

names(dur.ests.china) <- col.names
names(prev.ests.1) <- names(prev.ests.2) <- c(as.character(seq(y1,y2,1)),"Assumed duration of untreated disease","Duration of treated disease")
prev.ests.1.c <- prev.ests.1
prev.ests.2.c <- prev.ests.2

####################
##    CAMBODIA    ##
####################
country.tmp <- "Cambodia"
cam.perc.change <- 0.452
cam.prev1 <- dat.prev$`1st Prevalence`[dat.prev$Country==country.tmp]

cam.ann.rate <- getRateofChange(prev1=cam.prev1,prev2=817,year1=2002,year2=2010)

tmp.pop <- read.csv("Data/camb-pop.csv")


y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]

# run model for untreated TB assumed to be between 3 and 7 years
sim.dat.camb <- list()
dur.ests.camb <- data.frame(matrix(0,nr=length(c(3:7)),nc=4))
prev.ests.1 <- prev.ests.2 <- data.frame(matrix(0,nr=M*length(c(3:7)),nc=y2-y1+3))
for(j in 3:7){
  tmp.dat <- getDurationABC(country=country.tmp,pop.data=tmp.pop,duration.range=seq(1,6,0.01),unDxTBdur=j,M=M,ABC.cut.off=0.15)
  sim.dat.camb[[j-2]] <- tmp.dat
  
  # add estimates of duration
  dur.ests.camb[j-2,] <- c(tmp.dat$duration,tmp.dat$CIs,tmp.dat$Duration.u)
  
  # get prevalence estimates from each simulation
  tmp.prev <- getSummaryData(getDurationABC.output=tmp.dat)
  
  prev.ests.1[((j-3)*M+1):((j-3)*M+M),] <- cbind(tmp.prev$ann.prev.1,rep(j,M),rep(tmp.dat$duration,M))
  prev.ests.2[((j-3)*M+1):((j-3)*M+M),] <- cbind(tmp.prev$ann.prev.2,rep(j,M),rep(tmp.dat$duration,M))
  
}

names(dur.ests.camb) <- col.names
names(prev.ests.1) <- names(prev.ests.2) <- c(as.character(seq(y1,y2,1)),"Assumed duration of untreated disease","Duration of treated disease")
prev.ests.1.ca <- prev.ests.1
prev.ests.2.ca <- prev.ests.2

save.image()

##################
##   MYANMAR    ##
##################
country.tmp <- "Myanmar"
mya.perc.change <- 0.51
mya.prev1 <- dat.prev$`1st Prevalence`[dat.prev$Country==country.tmp]

mya.ann.rate <- getRateofChange(prev1=mya.prev1,prev2=415,year1=2009,year2=2017)

tmp.pop <- read.csv("Data/myanmar-pop.csv")

y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]

# run model for untreated TB assumed to be between 3 and 7 years
sim.dat.mya <- list()
dur.ests.mya <- data.frame(matrix(0,nr=length(c(3:7)),nc=4))
prev.ests.1 <- prev.ests.2 <- data.frame(matrix(0,nr=M*length(c(3:7)),nc=y2-y1+3))
for(j in 3:7){
  tmp.dat <- getDurationABC(country=country.tmp,pop.data=tmp.pop,duration.range=seq(1,6,0.01),unDxTBdur=j,M=M,ABC.cut.off=0.15)
  sim.dat.mya[[j-2]] <- tmp.dat
  
  # add estimates of duration
  dur.ests.mya[j-2,] <- c(tmp.dat$duration,tmp.dat$CIs,tmp.dat$Duration.u)
  
  # get prevalence estimates from each simulation
  tmp.prev <- getSummaryData(getDurationABC.output=tmp.dat)
  
  prev.ests.1[((j-3)*M+1):((j-3)*M+M),] <- cbind(tmp.prev$ann.prev.1,rep(j,M),rep(tmp.dat$duration,M))
  prev.ests.2[((j-3)*M+1):((j-3)*M+M),] <- cbind(tmp.prev$ann.prev.2,rep(j,M),rep(tmp.dat$duration,M))
  
}

names(dur.ests.mya) <- col.names
names(prev.ests.1) <- names(prev.ests.2) <- c(as.character(seq(y1,y2,1)),"Assumed duration of untreated disease","Duration of treated disease")
prev.ests.1.m <- prev.ests.1
prev.ests.2.m <- prev.ests.2

save(prev.ests.1.c,prev.ests.1.ca,prev.ests.1.i,prev.ests.1.m,prev.ests.1.p,prev.ests.1.v,file="Results/simPrevalences-1.RData")

save(prev.ests.2.c,prev.ests.2.ca,prev.ests.2.i,prev.ests.2.m,prev.ests.2.p,prev.ests.2.v,file="Results/simPrevalences-2.RData")

save(dur.ests.camb,dur.ests.china,dur.ests.indo,dur.ests.mya,dur.ests.phil,dur.ests.phil,dur.ests.viet,file="Results/durEstsABC.RData")

save(sim.dat.camb,sim.dat.china,sim.dat.indo,sim.dat.mya,sim.dat.phil,sim.dat.viet,file="allSimResultsABC.RData")

save.image()
