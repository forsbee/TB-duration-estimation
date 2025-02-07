###############################################################
## Program to run the country simulations using ABC approach ##
###############################################################

source("R/utils.R")
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
# prevalence survey estimates, with years and CIs-
## 10/1/24: use >15 year old data updated from Bob
dat.prev <- read_excel("Data/PrevalenceSurveyData.xlsx",sheet="Over15yrs")

# data on treatments 
dat.cure <- read_excel("Data/AnnualTreatments.xlsx",sheet="Sheet1")

# number of simulations to run
M <- 100000 # number of samples to generate
M2 <- 10000 # number to include in final empirical estimate of duration distribution

col.names <- c("Duration of treated disease","LCL","UCL","Assumed duration of untreated disease")

#############
## VIETNAM ##
#############
country.tmp <- "Vietnam"

viet.ann.rate <- 1-0.045

# get first prevalence survey value
## Bob: back-calculate from 2017-8 such that there is a 4.5% decline
prev1 <- dat.prev %>% filter(Country==country.tmp) %>% getPrev1(annRate = viet.ann.rate)
tmp.dat <- dat.prev[dat.prev$Country==country.tmp,]
error <- mean(c(tmp.dat$`2nd Survey Prevalence`-tmp.dat$`lower CI prev2`,
                tmp.dat$`upper CI prev2`-tmp.dat$`2nd Survey Prevalence`))

prev1.LCL <- prev1 - error
prev1.UCL <- prev1 + error
dat.prev[dat.prev$Country==country.tmp,]$`1st Prevalence` <- prev1
dat.prev[dat.prev$Country==country.tmp,]$`lower CI prev1` <- prev1.LCL
dat.prev[dat.prev$Country==country.tmp,]$`upper CI prev1` <- prev1.UCL
####

tmp.pop <- read.csv("Data/vietnam-pop.csv")
y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]

# run model for untreated TB assumed to be between 3 and 7 years
sim.dat.viet <- list()
dur.ests.viet<- data.frame(matrix(0,nr=length(c(3:7)),nc=4))
prev.ests <- data.frame(matrix(0,nr=M2*length(c(3:7)),nc=y2-y1+3))
for(j in 3:7){
  tmp.dat <- getDurationSIR(country=country.tmp,pop.data=tmp.pop,unDxTBdur=j,M=M,M2=M2)
  sim.dat.viet[[j-2]] <- tmp.dat
  
  # add estimates of duration
  dur.ests.viet[j-2,] <- c(tmp.dat$duration,tmp.dat$dur.CIs,tmp.dat$Duration.u)
  
  # get prevalence estimates from each simulation
  prev.ests[((j-3)*M2+1):((j-3)*M2+M2),] <- cbind(tmp.dat$resamp.Prevalences,rep(j,M2),rep(tmp.dat$duration,M2))
  
}

names(dur.ests.viet) <- col.names
names(prev.ests) <- c(as.character(seq(y1,y2,1)),"Assumed duration of untreated disease","Duration of treated disease")
prev.ests.v <- prev.ests

###############
## INDONESIA ##
###############
country.tmp <- "Indonesia"

# Per Bob, need to calculate first prevalence such that there is a 50% increase between surveys
# therefore prev1=prev2/1.5
tmp.dat <- dat.prev %>% filter(Country==country.tmp)
prev1 <- tmp.dat$`2nd Survey Prevalence`/1.5

error <- mean(c(tmp.dat$`2nd Survey Prevalence`-tmp.dat$`lower CI prev2`,
                tmp.dat$`upper CI prev2`-tmp.dat$`2nd Survey Prevalence`))
prev1.LCL <- prev1 - error
prev1.UCL <- prev1 + error

dat.prev[dat.prev$Country==country.tmp,]$`1st Prevalence` <- prev1
dat.prev[dat.prev$Country==country.tmp,]$`lower CI prev1` <- prev1.LCL
dat.prev[dat.prev$Country==country.tmp,]$`upper CI prev1` <- prev1.UCL

indo.ann.rate <- (tmp.dat$`2nd Survey Prevalence`/prev1)^(1/as.numeric(tmp.dat$`Time between 1st survey and 2nd survey (in years)`))-1

tmp.pop <- read.csv("Data/indo-pop.csv")

y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]

# run model for untreated TB assumed to be between 3 and 7 years
sim.dat.indo <- list()
dur.ests.indo<- data.frame(matrix(0,nr=length(c(3:7)),nc=4))
prev.ests <- data.frame(matrix(0,nr=M2*length(c(3:7)),nc=y2-y1+3))
for(j in 3:7){
  tmp.dat <- getDurationSIR(country=country.tmp,pop.data=tmp.pop,duration.range=seq(1,6,0.01),unDxTBdur=j,M=M,M2=M2)
  sim.dat.indo[[j-2]] <- tmp.dat
  
  # add estimates of duration
  dur.ests.indo[j-2,] <- c(tmp.dat$duration,tmp.dat$dur.CIs,tmp.dat$Duration.u)
  
  # get prevalence estimates from each simulation
  prev.ests[((j-3)*M2+1):((j-3)*M2+M2),] <- cbind(tmp.dat$resamp.Prevalences,rep(j,M2),rep(tmp.dat$duration,M2))
  
}

names(dur.ests.indo) <- col.names
names(prev.ests) <- c(as.character(seq(y1,y2,1)),"Assumed duration of untreated disease","Duration of treated disease")
prev.ests.i <- prev.ests

#################
## PHILIPPINES ##
#################
country.tmp <- "Philippines"

# Per Bob, need to calculate first prevalence such that there is a 10.5% increase between surveys
# therefore prev1=prev2/1.105
tmp.dat <- dat.prev %>% filter(Country==country.tmp)
prev1 <- tmp.dat$`2nd Survey Prevalence`/1.105

error <- mean(c(tmp.dat$`2nd Survey Prevalence`-tmp.dat$`lower CI prev2`,
                tmp.dat$`upper CI prev2`-tmp.dat$`2nd Survey Prevalence`))
prev1.LCL <- prev1 - error
prev1.UCL <- prev1 + error

dat.prev[dat.prev$Country==country.tmp,]$`1st Prevalence` <- prev1
dat.prev[dat.prev$Country==country.tmp,]$`lower CI prev1` <- prev1.LCL
dat.prev[dat.prev$Country==country.tmp,]$`upper CI prev1` <- prev1.UCL

phil.ann.rate <- (tmp.dat$`2nd Survey Prevalence`/prev1)^(1/as.numeric(tmp.dat$`Time between 1st survey and 2nd survey (in years)`))-1
####

tmp.pop <- read.csv("Data/phil-pop.csv")

y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]

# run model for untreated TB assumed to be between 3 and 7 years
sim.dat.phil <- list()
dur.ests.phil <- data.frame(matrix(0,nr=length(c(3:7)),nc=4))
prev.ests <- data.frame(matrix(0,nr=M2*length(c(3:7)),nc=y2-y1+3))
for(j in 3:7){
  tmp.dat <- getDurationSIR(country=country.tmp,pop.data=tmp.pop,duration.range=seq(1,6,0.01),unDxTBdur=j,M=M,M2=M2)
  sim.dat.phil[[j-2]] <- tmp.dat
  
  # add estimates of duration
  dur.ests.phil[j-2,] <- c(tmp.dat$duration,tmp.dat$dur.CIs,tmp.dat$Duration.u)
  
  # get prevalence estimates from each simulation
  prev.ests[((j-3)*M2+1):((j-3)*M2+M2),] <- cbind(tmp.dat$resamp.Prevalences,rep(j,M2),rep(tmp.dat$duration,M2))
  
}

names(dur.ests.phil) <- col.names
names(prev.ests) <- c(as.character(seq(y1,y2,1)),"Assumed duration of untreated disease","Duration of treated disease")
prev.ests.p <- prev.ests

#################
##    CHINA    ##
#################
country.tmp <- "China"
china.ann.rate <- dat.prev %>% filter(Country==country.tmp) %>% getRateofChange2()

tmp.pop <- read.csv("Data/china-pop.csv")

y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]

# run model for untreated TB assumed to be between 3 and 7 years
sim.dat.china <- list()
dur.ests.china <- data.frame(matrix(0,nr=length(c(3:7)),nc=4))
prev.ests <- data.frame(matrix(0,nr=M2*length(c(3:7)),nc=y2-y1+3))
for(j in 3:7){
  tmp.dat <- getDurationSIR(country=country.tmp,pop.data=tmp.pop,duration.range=seq(1,6,0.01),unDxTBdur=j,M=M,M2=M2)
  sim.dat.china[[j-2]] <- tmp.dat
  
  # add estimates of duration
  dur.ests.china[j-2,] <- c(tmp.dat$duration,tmp.dat$dur.CIs,tmp.dat$Duration.u)
  
  # get prevalence estimates from each simulation
  prev.ests[((j-3)*M2+1):((j-3)*M2+M2),] <- cbind(tmp.dat$resamp.Prevalences,rep(j,M2),rep(tmp.dat$duration,M2))
  
}

names(dur.ests.china) <- col.names
names(prev.ests) <- c(as.character(seq(y1,y2,1)),"Assumed duration of untreated disease","Duration of treated disease")
prev.ests.c <- prev.ests

####################
##    CAMBODIA    ##
####################
country.tmp <- "Cambodia"
cam.ann.rate <- dat.prev %>% filter(Country==country.tmp) %>% getRateofChange2()

tmp.pop <- read.csv("Data/camb-pop.csv")

y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]

# run model for untreated TB assumed to be between 3 and 7 years
sim.dat.camb <- list()
dur.ests.camb <- data.frame(matrix(0,nr=length(c(3:7)),nc=4))
prev.ests <- data.frame(matrix(0,nr=M2*length(c(3:7)),nc=y2-y1+3))
for(j in 3:7){
  tmp.dat <- getDurationSIR(country=country.tmp,pop.data=tmp.pop,duration.range=seq(1,6,0.01),unDxTBdur=j,M=M,M2=M2)
  sim.dat.camb[[j-2]] <- tmp.dat
  
  # add estimates of duration
  dur.ests.camb[j-2,] <- c(tmp.dat$duration,tmp.dat$dur.CIs,tmp.dat$Duration.u)
  
  # get prevalence estimates from each simulation
  prev.ests[((j-3)*M2+1):((j-3)*M2+M2),] <- cbind(tmp.dat$resamp.Prevalences,rep(j,M2),rep(tmp.dat$duration,M2))
  
}

names(dur.ests.camb) <- col.names
names(prev.ests) <- c(as.character(seq(y1,y2,1)),"Assumed duration of untreated disease","Duration of treated disease")
prev.ests.ca <- prev.ests

##################
##   MYANMAR    ##
##################
country.tmp <- "Myanmar"

# total change between first and second prevalence is 51%; therefore first prevalence should be prev2/1.51
tmp.dat <- dat.prev %>% filter(Country==country.tmp)
prev1 <- tmp.dat$`2nd Survey Prevalence`/0.49

error <- mean(c(tmp.dat$`2nd Survey Prevalence`-tmp.dat$`lower CI prev2`,
                tmp.dat$`upper CI prev2`-tmp.dat$`2nd Survey Prevalence`))
prev1.LCL <- prev1 - error
prev1.UCL <- prev1 + error

dat.prev[dat.prev$Country==country.tmp,]$`1st Prevalence` <- prev1
dat.prev[dat.prev$Country==country.tmp,]$`lower CI prev1` <- prev1.LCL
dat.prev[dat.prev$Country==country.tmp,]$`upper CI prev1` <- prev1.UCL

mya.ann.rate <- 1-(tmp.dat$`2nd Survey Prevalence`/prev1)^(1/as.numeric(tmp.dat$`Time between 1st survey and 2nd survey (in years)`))

#####

tmp.pop <- read.csv("Data/myanmar-pop.csv")

y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]

# run model for untreated TB assumed to be between 3 and 7 years
sim.dat.mya <- list()
dur.ests.mya <- data.frame(matrix(0,nr=length(c(3:7)),nc=4))
prev.ests <- data.frame(matrix(0,nr=M2*length(c(3:7)),nc=y2-y1+3))
for(j in 3:7){
  tmp.dat <- getDurationSIR(country=country.tmp,pop.data=tmp.pop,duration.range=seq(1,6,0.01),unDxTBdur=j,M=M,M2=M2)
  sim.dat.mya[[j-2]] <- tmp.dat
  
  # add estimates of duration
  dur.ests.mya[j-2,] <- c(tmp.dat$duration,tmp.dat$dur.CIs,tmp.dat$Duration.u)
  
  # get prevalence estimates from each simulation
  prev.ests[((j-3)*M2+1):((j-3)*M2+M2),] <- cbind(tmp.dat$resamp.Prevalences,rep(j,M2),rep(tmp.dat$duration,M2))
  
}

names(dur.ests.mya) <- col.names
names(prev.ests) <- c(as.character(seq(y1,y2,1)),"Assumed duration of untreated disease","Duration of treated disease")
prev.ests.m <- prev.ests

save(prev.ests.c,prev.ests.ca,prev.ests.i,prev.ests.m,prev.ests.p,prev.ests.v,file="Results/simPrevalences_SIR.RData")

save(dur.ests.camb,dur.ests.china,dur.ests.indo,dur.ests.mya,dur.ests.phil,dur.ests.phil,dur.ests.viet,file="Results/durEstsSIR.RData")

save(sim.dat.camb,sim.dat.china,sim.dat.indo,sim.dat.mya,sim.dat.phil,sim.dat.viet,file="Results/allSimResultsSIR.RData")

save.image()
