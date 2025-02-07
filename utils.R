# FUNCTIONS FOR PREVALENCE CALCULATIONS #
#########################################
#function to get rate of change assuming multiplicative rate of change
getRateofChange <- function(prev1,prev2,year1,year2){
  time.period <- year2-year1
  r <- (prev2/prev1)^(1/time.period)
  
  return(r)
}

#function to get rate of change assuming multiplicative rate of change
getRateofChange2 <- function(dat.prev.country){

   time.period <- dat.prev.country$`Second prevalence survey year` - dat.prev.country$`First prevalence survey year`
  r <- (dat.prev.country$`2nd Survey Prevalence`/dat.prev.country$`1st Prevalence`)^(1/time.period)
  
  return(r)
}

#function to get rate of change assuming multiplicative rate of change
getPrev1 <- function(dat.prev.country,annRate){
  
  time.period <- dat.prev.country$`Second prevalence survey year` - dat.prev.country$`First prevalence survey year`
  prev1 <- dat.prev.country$`2nd Survey Prevalence`/(annRate^time.period)
  
  return(prev1)
}


#########################################################################################################################
# MAIN MODEL FUNCTION
######
# function that takes a ratio between incidence and prevalence as input and calculates annual prevalence from this
## used to search over this ratio to find the optimal value to arrive at second prevalence value that is closest to actual value
calcAnnPrev <- function(duration,prev1,cures,years,unDxTBdur=5){
  prevs <- NULL
  prevs[1] <- prev1
  
  for(i in 1:(length(years)-1)){
    prevs[i+1] <- prevs[i] - cures[i] - (1/unDxTBdur)*(prevs[i] - cures[i]) + (1/duration)*prevs[i]
  }
  
  return(prevs)
}

#########################################################################################################################
# SAMPLING IMPORTANCE RESAMPLING ESTIMATION ALGORITHM 
#####
# use sample importance resampling approach to get posterior estimate of duration of untreated disease
# assumes that there is an object data.prev and dat.cure which has prevalence survey data (in runSIRsims.R)
getDurationSIR <- function(country,pop.data,duration.range=seq(1,6,0.01),unDxTBdur=5,M=100000,M2=10000){
  set.seed(21)
  
  # get distribution of P_1, and P_2 (first and second prevalences)-we assume they are normally distributed
  dat.prev.tmp <- subset(dat.prev,dat.prev$Country==country)
  sd.1.l <- (dat.prev.tmp$`1st Prevalence`-dat.prev.tmp$`lower CI prev1`)/1.96
  sd.1.u <- abs(dat.prev.tmp$`1st Prevalence`-dat.prev.tmp$`upper CI prev1`)/1.96
  sd.1 <- mean(c(sd.1.l,sd.1.u))
  
  sd.2.l <- (dat.prev.tmp$`2nd Survey Prevalence`-dat.prev.tmp$`lower CI prev2`)/1.96
  sd.2.u <- abs(dat.prev.tmp$`2nd Survey Prevalence`-dat.prev.tmp$`upper CI prev2`)/1.96
  sd.2 <- mean(c(sd.2.l,sd.2.u))
  
  # get years of surveys
  years.tmp <- c(dat.prev.tmp$`First prevalence survey year`:dat.prev.tmp$`Second prevalence survey year`)
  y1 <- years.tmp[1]
  y2 <- years.tmp[length(years.tmp)]
  
  # get first and second population size estimates, in  thousands
  pop.yrs <- pop.data$Total.Population[which(pop.data$label==y1):which(pop.data$label==y2)]
  pop1 <- pop.yrs[1]
  pop2 <- pop.yrs[length(years.tmp)]
  
  # get number of cases per 1000
  prev1.n <- (dat.prev.tmp$`1st Prevalence`/1000)*(pop1/100000)
  
  # extract the number of cures
  cures.tmp <- dat.cure$cured[dat.cure$Country==country] # number of cures for the country (in thousands)

  # now run the simulation
  # generate M prev1 values from normal distribution (per 100k)
  prev.1 <- rnorm(M,mean=dat.prev.tmp$`1st Prevalence`,sd=sd.1)
  
  # convert these to raw number of cases to run the model to get prev2
  prev.1.n <- (prev.1/100000)*(pop1/1000)
  
  # sample M disease durations
  dur.tmp <- runif(M,1,6)
  
  # calculate resulting number of cases at the second prevalence survey time
  prev.n.2 <- NULL
  for(i in 1:M){
    prev.n.2[i] <- calcAnnPrev(duration=dur.tmp[i],prev1=prev.1.n[i],cures=cures.tmp,years=years.tmp,unDxTBdur)[length(years.tmp)]
  }
  
  # need to convert back to population adjusted to match distn above
  prev.2 <- (1000*prev.n.2/pop2)*100000

  # calculate the probability of observing prev.1 and prev.2 using normal distribution
  p.1 <- dnorm((prev.1-dat.prev.tmp$`1st Prevalence`)/sd.1)
  p.2 <- dnorm((prev.2-dat.prev.tmp$`2nd Survey Prevalence`)/sd.2)
  
  # remove those inds for which p.2<=0.0001 to ensure non-zero values
  inds.keep <- which(p.2 >= 0.0001)
  dur.est.tmp <- dur.tmp[inds.keep] # keep duration values that allow p.2>= 0.01)
  p1.keep <- prev.1[inds.keep] # prevalence 1 values that keep prev 2 above threshold
  prev.1.n.keep <- prev.1.n[inds.keep] # prevalence 1 in terms of number of cases
  
  # calculate a weight for this as product of p.1 to p.2
  weights.tmp <- p.1[inds.keep]*p.2[inds.keep] 

  # normalize weights
  weights.norm <- weights.tmp/sum(weights.tmp)
  
  # sample with replacement from durations according to normalized weights
  inds.dist <- sample(c(1:length(weights.norm)),size=M2,replace=T,prob=weights.norm)
  duration.dist <- dur.est.tmp[inds.dist]
  p1.dist <- p1.keep[inds.dist]
  prev.1.n.dist <- prev.1.n.keep[inds.dist]

  # now get estimate of duration and empirical CIs
  dur.est <- mean(duration.dist)
  dur.CIs <- quantile(duration.dist,c(0.025,0.975))
  
  # get corresponding simulated annual prevalence values (to be used to get CIs for prevalence)
  sim.n.prevs <- matrix(0,nr=M2,nc=length(years.tmp))
  for(i in 1:M2){
    sim.n.prevs[i,] <- calcAnnPrev(duration=duration.dist[i],prev1=prev.1.n.dist[i],
                                   cures=cures.tmp,years=years.tmp,unDxTBdur=unDxTBdur)
  }
  # convert to population-adjusted values
  sim.prevs <- t((10^8)*(t(sim.n.prevs)/pop.yrs))
  
  # get 95% CIs of annual prevalences from these
  prev.CIs.sim <- apply(sim.prevs,2,quantile,c(0.025,0.975))
  prev.mean.sim <- apply(sim.prevs,2,mean)
  
  return(list(Country=country,duration=dur.est,dur.CIs=dur.CIs,
              Years=seq(y1,y2,1),
              resamp.Prevalences=sim.prevs,resamp.prev.est=prev.mean.sim,resamp.prev.CIs=prev.CIs.sim,
              Treated=cures.tmp,sim.Durations=duration.dist,Duration.u=unDxTBdur))
}

#########################################################################################################
# function to convert prevalence data from CalcAnnPrev to long format
getPrevLong <- function(prev.dat.wide,y1,y2,pop.dat,cure.dat,country.name,simDataABC){
  
  names(pop.dat) <- c("Year","Popn")
  
  prev.dat.wide <- data.frame(prev.dat.wide)
  names(prev.dat.wide) <- c("3","4","5","6","7")
  years <- data.frame(as.character(c(y1:y2)))
  names(years) <- "Year"
  prev.dat.wide <- bind_cols(years,prev.dat.wide)
  
  # transform to long format
  long.prev.dat <- prev.dat.wide %>%
    pivot_longer(!Year, names_to = "Untreated duration", values_to = "Prevalence")
  
  long.prev.dat.pop <- merge(long.prev.dat, pop.dat, by = "Year") %>% 
    mutate(Prev100k = (10^8)*Prevalence/Popn)
  
  Country <- data.frame(rep(country.name,dim(long.prev.dat.pop)[1]))
  names(Country) <- "Country"
  
  long.prev.dat.pop <- bind_cols(Country,long.prev.dat.pop)
  
  # add in proportion of prevalent cases treated
  # first add in cure data to the data frame
  cures.tmp <- data.frame(Year=c((y1+1):y2),Cures=cure.dat)
  long.dat <- merge(long.prev.dat.pop,cures.tmp,all.x=T) 
  
  # now calculate proportion treated each year
  long.dat$`Prop treated` <- long.dat$Cures/long.dat$Prevalence
  
  # add in CIs
  l <- length(c(y1:y2))
  CIs <- data.frame(`Untreated duration`=rep(c(3:7),each=l),Year=rep(c(y1:y2),times=5),
                    LCL=rep(NA,length(long.dat[,1])),UCL=rep(NA,length(long.dat[,1])),
                    Prev100k.sim=rep(NA,length(long.dat[,1])),
                    PrevDelta=rep(NA,length(long.dat[,1])))
  for(i in 1:length(simDataABC)){
    CIs$LCL[(l*(i-1)+1):(l*i)] <- simDataABC[[i]]$resamp.prev.CIs[1,]
    CIs$UCL[(l*(i-1)+1):(l*i)] <- simDataABC[[i]]$resamp.prev.CIs[2,]
    CIs$Prev100k.sim[(l*(i-1)+1):(l*i)] <- simDataABC[[i]]$resamp.prev.est
    CIs$PrevDelta[(l*(i-1)+1):(l*i)] <- c(NA,diff(simDataABC[[i]]$resamp.prev.est))
  }
  CIs <- CIs[order(CIs$Untreated.duration),]
  names(CIs) <- c("Untreated duration","Year","LCL","UCL","Prev100k.sim","PrevDelta")
  
  # merge this into rest of the data
  long.dat.all <- merge(long.dat,CIs)
  
  return(long.dat.all)
}

########################################################################################################################
# deprecated functions below here
#########################################################################################################################

# function to account for uncertainty in first and second prevalence survey values
## uses calcAnnPrev() and samples first prevalence from distn to capture uncertainty in that value
## then it calculates prevalence up to the second prevalence survey value for a range of duration values
## duration values are selected if they produce p_1 and p_T values that meet criteria for acceptance
## based on assuming they are normally distributed
getDurationABC_mv <- function(country,pop.data,duration.range=seq(1,6,0.01),unDxTBdur=5,M=1000,ABC.cut.off=0.5){
  # This approach uses weights to attempt to control variance through the vectors of probabilities
  set.seed(21)
  
  # get distribution of P_1, and P_2 (first and second prevalences)-we assume they are normally distributed
  dat.prev.tmp <- subset(dat.prev,dat.prev$Country==country)
  sd.1.l <- (dat.prev.tmp$`1st Prevalence`-dat.prev.tmp$`lower CI prev1`)/1.96
  sd.1.u <- abs(dat.prev.tmp$`1st Prevalence`-dat.prev.tmp$`upper CI prev1`)/1.96
  sd.1 <- mean(c(sd.1.l,sd.1.u))
  
  sd.2.l <- (dat.prev.tmp$`2nd Survey Prevalence`-dat.prev.tmp$`lower CI prev2`)/1.96
  sd.2.u <- abs(dat.prev.tmp$`2nd Survey Prevalence`-dat.prev.tmp$`upper CI prev2`)/1.96
  sd.2 <- mean(c(sd.2.l,sd.2.u))
  
  # get years of surveys
  years.tmp <- c(dat.prev.tmp$`First prevalence survey year`:dat.prev.tmp$`Second prevalence survey year`)
  y1 <- years.tmp[1]
  y2 <- years.tmp[length(years.tmp)]
  
  # get first and second population size estimates, in  thousands
  pop.yrs <- pop.data$Total.Population[which(pop.data$label==y1):which(pop.data$label==y2)]
  pop1 <- pop.yrs[1]
  pop2 <- pop.yrs[length(years.tmp)]
  
  # get number of cases per 1000
  prev1.n <- (dat.prev.tmp$`1st Prevalence`/1000)*(pop1/100000)
  
  # now run the simulation
  cures.tmp <- dat.cure$cured[dat.cure$Country==country] # number of cures for the country (in thousands)
  z.scores <- data.frame(matrix(nr=0,nc=length(duration.range))) #space to store z scores used to select durations
  dur.est.accept <- NULL #vector to store duration selected for each simulation
  prev.n.sim.vecs <- prev.sim.vecs <- data.frame(matrix(ncol=length(years.tmp),nrow=0))
  counter <- 1 # counts number of accepted samples
  counter.2 <- 1 # counts number of loops that need to run to get M accepted values
  while(counter<=M){
    # generate a first prevalence survey value (per 100k)
    prev.1 <- rnorm(1,mean=dat.prev.tmp$`1st Prevalence`,sd=sd.1)
    
    # convert this to the raw number of cases for running the model
    prev.1.n <- (prev.1/100000)*(pop1/1000) 
    
    # get list with each item being a vector of prevalences over the years between surveys for each proposed duration value
    prevs.tmp.n <- lapply(duration.range,calcAnnPrev,prev1=prev.1.n,cures=cures.tmp,years=years.tmp,unDxTBdur=unDxTBdur)
    # convert this list to a matrix with as many rows as the length of duration.range
    prevs.tmp.n.2 <- matrix(unlist(prevs.tmp.n),nr=length(duration.range),byrow=T)
    
    # Get last estimated prevalence for each of these
    prev.2.n <- prevs.tmp.n.2[,dim(prevs.tmp.n.2)[2]]
    
    # need to get these prevalences back to being population adjusted to match distn above
    prev.2.tmp <- (1000*prev.2.n/pop2)*100000
    
    # calculate the probability of observing prev.1 and prev.2 using mv normal distribution (Assuming independence between p1 and p2)
    p.1 <- 2*(1-pnorm(abs((prev.1-dat.prev.tmp$`1st Prevalence`)/sd.1)))
    p.2 <- 2*(1-pnorm(abs((prev.2.tmp-dat.prev.tmp$`2nd Survey Prevalence`)/sd.2)))
    p.vec <- p.1*p.2
    
    # keep indices of durations with probabilities that are greater than than ABC.cut.off
    inds.accept <- which(p.vec>ABC.cut.off)
    
    if(length(inds.accept)>0){
      # keep prev.tmp.n that correspond to durations that we keep
      prev.n.sim.vecs[counter:(counter+length(inds.accept)-1),] <- data.frame(matrix(unlist(prevs.tmp.n[inds.accept]),nc=length(years.tmp),byrow=T))
      
      # store the duration estimate(s) that were selected
      dur.est.accept[counter:(counter+length(inds.accept)-1)] <- duration.range[inds.accept]
      
      counter <- counter+length(inds.accept)
    }
    
    counter.2 <- counter.2+1
    
  }
  
  # get prevalence that is population adjusted
  prev.sim.vecs <- t((10^8)*(t(prev.n.sim.vecs)/pop.yrs))
  
  # now get estimate of duration and empirical CIs
  dur.est <- mean(dur.est.accept)
  dur.CIs <- quantile(dur.est.accept,c(0.025,0.975))
  
  # get 95% CIs and mean for prevalence estimates as counts
  prev.n.mean <- apply(prev.n.sim.vecs,2,mean)
  prev.n.CIs <- apply(prev.n.sim.vecs,2,quantile,c(0.025,0.975))
  
  # get 95% CIs and mean for population-adjusted prevalence estimates
  prev.mean <- apply(prev.sim.vecs,2,mean)
  prev.CIs <- apply(prev.sim.vecs,2,quantile,c(0.025,0.975))
  
  # get 95% CIs using more accurate resampling approach
  sim.n.prevs <- matrix(0,nr=M,nc=length(years.tmp))
  # randomly select M first prevalences and durations to simulate from
  inds <- sample(c(1:M),replace=T)
  
  for(i in 1:M){
    sim.n.prevs[i,] <- calcAnnPrev(duration=dur.est.accept[inds[i]],prev1=prev.n.sim.vecs[inds[i],1],
                                   cures=cures.tmp,years=years.tmp,unDxTBdur=unDxTBdur)
  }
  # convert to population-adjusted values
  sim.prevs <- t((10^8)*(t(sim.n.prevs)/pop.yrs))
  
  # get 95% CIs from these
  prev.CIs.sim <- apply(sim.prevs,2,quantile,c(0.025,0.975))
  prev.mean.sim <- apply(sim.prevs,2,mean)
  
  return(list(Country=country,duration=dur.est,dur.CIs=dur.CIs,Years=seq(y1,y2,1),
              sim.Prevalences=prev.sim.vecs,prev.est=prev.mean,prev.CIs=prev.CIs,
              sim.n.Prevalences=prev.n.sim.vecs,prev.n.est=prev.n.mean,prev.n.CIs=prev.n.CIs,
              resamp.Prevalences=sim.prevs,resamp.prev.est=prev.mean.sim,resamp.prev.CIs=prev.CIs.sim,
              Treated=cures.tmp,sim.Durations=dur.est.accept,Duration.u=unDxTBdur,z.scores=z.scores,counter=counter.2))
}

# function to account for uncertainty in first and second prevalence survey values
## uses calcAnnPrev() and samples first prevalence from distn to capture uncertainty in that value
## then it assigns a weight to each potential duration based on the pdf of the second prevalence value
# duration values are accepted or rejected based on that weight
getDurationABC <- function(country,pop.data,duration.range=seq(1,6,0.01),unDxTBdur=5,M=1000,ABC.cut.off=qnorm(0.525)){
  set.seed(21)

  # get distribution of P_1, and P_2 (first and second prevalences)-we assume they are normally distributed
  dat.prev.tmp <- subset(dat.prev,dat.prev$Country==country)
  sd.1.l <- (dat.prev.tmp$`1st Prevalence`-dat.prev.tmp$`lower CI prev1`)/1.96
  sd.1.u <- abs(dat.prev.tmp$`1st Prevalence`-dat.prev.tmp$`upper CI prev1`)/1.96
  sd.1 <- mean(c(sd.1.l,sd.1.u))
  
  sd.2.l <- (dat.prev.tmp$`2nd Survey Prevalence`-dat.prev.tmp$`lower CI prev2`)/1.96
  sd.2.u <- abs(dat.prev.tmp$`2nd Survey Prevalence`-dat.prev.tmp$`upper CI prev2`)/1.96
  sd.2 <- mean(c(sd.2.l,sd.2.u))

  # get years of surveys
  years.tmp <- c(dat.prev.tmp$`First prevalence survey year`:dat.prev.tmp$`Second prevalence survey year`)
  y1 <- years.tmp[1]
  y2 <- years.tmp[length(years.tmp)]
  
  # get first and second population size estimates, in  thousands
  pop.yrs <- pop.data$Total.Population[which(pop.data$label==y1):which(pop.data$label==y2)]
  pop1 <- pop.yrs[1]
  pop2 <- pop.yrs[length(years.tmp)]
  
  # get number of cases per 1000
  prev1.n <- (dat.prev.tmp$`1st Prevalence`/1000)*(pop1/100000)
  
  # now run the simulation
  cures.tmp <- dat.cure$cured[dat.cure$Country==country] # number of cures for the country (in thousands)
  z.scores <- data.frame(matrix(nr=0,nc=length(duration.range))) #space to store z scores used to select durations
  dur.est.accept <- NULL #vector to store duration selected for each simulation
  prev.n.sim.vecs <- prev.sim.vecs <- data.frame(matrix(ncol=length(years.tmp),nrow=0))
  counter <- 1 # counts number of accepted samples
  counter.2 <- 1 # counts number of loops that need to run to get M accepted values
  while(counter<=M){
    # generate a first prevalence survey value (per 100k)
    prev.1 <- rnorm(1,mean=dat.prev.tmp$`1st Prevalence`,sd=sd.1)
    # convert this to the raw number of cases for running the model
    prev.1.n <- (prev.1/100000)*(pop1/1000) 
    
    # get list with each item being a vector of prevalences over the years between surveys for each proposed duration value
    prevs.tmp.n <- lapply(duration.range,calcAnnPrev,prev1=prev.1.n,cures=cures.tmp,years=years.tmp,unDxTBdur=unDxTBdur)
    # convert this list to a matrix with as many rows as the length of duration.range
    prevs.tmp.n.2 <- matrix(unlist(prevs.tmp.n),nr=length(duration.range),byrow=T)
    
    # Get last estimated prevalence for each of these
    prev.2.n <- prevs.tmp.n.2[,dim(prevs.tmp.n.2)[2]]
    
    # need to get these prevalences back to being population adjusted to match distn above
    prev.2.tmp <- (1000*prev.2.n/pop2)*100000
    
    # get z scores of prev.2.tmp values assuming second prevalence survey value follows normal distn
    z.scores[counter.2,] <- (prev.2.tmp-dat.prev.tmp$`2nd Survey Prevalence`)/sd.2
    
    # keep indices of durations with z.scores that have probability less than ABC.cut.off
    inds.accept <- getZ(z.scores[counter.2,],ABC.cut.off)
    
    if(length(inds.accept)>0){
      # keep prev.tmp.n that correspond to durations that we keep
      prev.n.sim.vecs[counter:(counter+length(inds.accept)-1),] <- data.frame(matrix(unlist(prevs.tmp.n[inds.accept]),nc=length(years.tmp),byrow=T))
      
      # store the duration estimate(s) that were selected
      dur.est.accept[counter:(counter+length(inds.accept)-1)] <- duration.range[inds.accept]
      
      counter <- counter+length(inds.accept)
    }
    
    counter.2 <- counter.2+1

  }
  
  # get prevalence that is population adjusted
  prev.sim.vecs <- t((10^8)*(t(prev.n.sim.vecs)/pop.yrs))
  
  # now get estimate of duration and empirical CIs
  dur.est <- mean(dur.est.accept)
  dur.CIs <- quantile(dur.est.accept,c(0.025,0.975))
  
  # get 95% CIs and mean for prevalence estimates as counts
  prev.n.mean <- apply(prev.n.sim.vecs,2,mean)
  prev.n.CIs <- apply(prev.n.sim.vecs,2,quantile,c(0.025,0.975))

  # get 95% CIs and mean for population-adjusted prevalence estimates
  prev.mean <- apply(prev.sim.vecs,2,mean)
  prev.CIs <- apply(prev.sim.vecs,2,quantile,c(0.025,0.975))
  
  # get 95% CIs using more accurate resampling approach
  sim.n.prevs <- matrix(0,nr=M,nc=length(years.tmp))
  # randomly select M first prevalences and durations to simulate from
  inds <- sample(c(1:M),replace=T)
  
  for(i in 1:M){
    sim.n.prevs[i,] <- calcAnnPrev(duration=dur.est.accept[inds[i]],prev1=prev.n.sim.vecs[inds[i],1],
                                  cures=cures.tmp,years=years.tmp,unDxTBdur=unDxTBdur)
  }
  # convert to population-adjusted values
  sim.prevs <- t((10^8)*(t(sim.n.prevs)/pop.yrs))
  
  # get 95% CIs from these
  prev.CIs.sim <- apply(sim.prevs,2,quantile,c(0.025,0.975))
  prev.mean.sim <- apply(sim.prevs,2,mean)
  
  return(list(Country=country,duration=dur.est,dur.CIs=dur.CIs,Years=seq(y1,y2,1),
              sim.Prevalences=prev.sim.vecs,prev.est=prev.mean,prev.CIs=prev.CIs,
              sim.n.Prevalences=prev.n.sim.vecs,prev.n.est=prev.n.mean,prev.n.CIs=prev.n.CIs,
              resamp.Prevalences=sim.prevs,resamp.prev.est=prev.mean.sim,resamp.prev.CIs=prev.CIs.sim,
              Treated=cures.tmp,sim.Durations=dur.est.accept,Duration.u=unDxTBdur,z.scores=z.scores,counter=counter.2))
}

#########################################################################################################################

# helper function to getDurationABC() to select positions in a vector where the value is less than a cut.off
getZ <- function(vec,cut.off){
  return(which(abs(vec) < cut.off))
}

#########################################################################################################################
getIncidenceEsts <- function(getDurationABC.output){
  # function to take output from getDurationABC() and calculate the associated incidence estimates
  
  # extract the accepted durations and prevalence vectors associated with this
  sim.Durations <- getDurationABC.output$sim.Durations
  prev.vecs <- getDurationABC.output$resamp.Prevalences
  
  # get incidence at all time points as: I=Prev/duration
  incidence.vecs <- prev.vecs/sim.Durations
  
  # get mean and CIs
  incidence.LCL <- apply(incidence.vecs,2,quantile,0.025)
  incidence.UCL <- apply(incidence.vecs,2,quantile,0.975)
  
  incidence.est <- apply(incidence.vecs,2,mean)
  
  return(list(incidence.est,incidence.LCL,incidence.UCL))
}

