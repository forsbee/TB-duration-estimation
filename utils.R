
#function to get rate of change assuming multiplicative rate of change
getRateofChange <- function(prev1,prev2,year1,year2){
  time.period <- year2-year1
  r <- (prev2/prev1)^(1/time.period)
  
  return(r)
}

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

# helper function to select postions in a vector where the value is less than a cut.off
getZ <- function(vec,cut.off){
  return(which(abs(vec) < cut.off))
}

# function to account for uncertainty in first and second prevalence survey values
## uses calcAnnPrev() and samples first prevalence from distn to capture uncertainty in that value
## then it assigns a weight to each potential duration based on the pdf of the second prevalence value
# duration values are accepted or rejected based on that weight
getDurationABC <- function(country,pop.data,duration.range=seq(1,6,0.01),unDxTBdur=5,M=1000,ABC.cut.off=0.15){
  # get distribution of P_1, and P_2 (first and second prevalences)-we assume they are normally distributed
  set.seed(21)
  dat.prev.tmp <- subset(dat.prev,dat.prev$Country==country)
  sd.1.l <- (dat.prev.tmp$`1st Prevalence`-dat.prev.tmp$`lower CI prev1`)/1.96
  sd.1.u <- abs(dat.prev.tmp$`1st Prevalence`-dat.prev.tmp$`upper CI prev1`)/1.96
  sd.1 <- mean(c(sd.1.l,sd.1.u))
  
  sd.2.l <- (dat.prev.tmp$`2nd Survey Prevalence`-dat.prev.tmp$`lower CI prev2`)/1.96
  sd.2.u <- abs(dat.prev.tmp$`2nd Survey Prevalence`-dat.prev.tmp$`upper CI prev2`)/1.96
  sd.2 <- mean(c(sd.2.l,sd.2.u))

  # get years of surveys
  y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country]
  y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country]
  
  # get first and second population size estimates, in  thousands
  pop1 <- pop.data$Total.Population[pop.data$label==y1]
  pop2 <- pop.data$Total.Population[pop.data$label==y2]
  
  # get number of cases per 1000
  prev1.n <- (dat.prev.tmp$`1st Prevalence`/1000)*(pop1/100000)
  
  # now run the simulation
  cures.tmp <- dat.cure$cured[dat.cure$Country==country] # number of cures for the country (in thousands)
  years.tmp <- dat.cure$Year[dat.cure$Country==country] # years of cure data
  z.scores <- data.frame(matrix(nr=0,nc=length(duration.range))) #space to store z scores used to select durations
  dur.est.accept <- NULL #vector to store duration selected for each simulation
  prev.sim.vecs <- data.frame(matrix(ncol=length(years.tmp),nrow=0))
  counter <- 1 # counts number of accepted samples
  counter.2 <- 0 # counts number of loops that need to run to get M accepted values
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
      prev.sim.vecs[counter:(counter+length(inds.accept)-1),] <- data.frame(matrix(unlist(prevs.tmp.n[inds.accept]),nc=length(years.tmp),byrow=T))
      
      # store the duration estimate(s) that were selected
      dur.est.accept[counter:(counter+length(inds.accept)-1)] <- duration.range[inds.accept]
      
      counter <- counter+length(inds.accept)
    }
    
    counter.2 <- counter.2+1

  }
  
  # now get estimate of duration and empirical CIs
  dur.est <- mean(dur.est.accept)
  CIs <- quantile(dur.est.accept,c(0.025,0.975))
    
  return(list(Country=country,duration=dur.est,CIs=CIs,Years=seq(y1,y2,1),sim.Prevalences=prev.sim.vecs,
              Treated=cures.tmp,sim.Durations=dur.est.accept,Duration.u=unDxTBdur,z.scores=z.scores))
}

getIncidenceEsts <- function(getDurationABC.output){
  # function to take output from getDurationABC() and calculate the associated incidence estimates
  
  # extract the accepted durations and prevalence vectors associated with this
  sim.Durations <- getDurationABC.output$sim.Durations
  prev.vecs <- getDurationABC.output$sim.Prevalences
  
  # get incidence at all time points as: I=Prev/duration
  incidence.vecs <- prev.vecs/sim.Durations
  
  # get mean and CIs
  incidence.LCL <- apply(incidence.vecs,2,quantile,0.025)
  incidence.UCL <- apply(incidence.vecs,2,quantile,0.975)
  
  incidence.est <- apply(incidence.vecs,2,mean)
  
  return(list(incidence.est,incidence.LCL,incidence.UCL))
}

# function to make nice outputs for prevalence and incidence
getSummaryData <- function(getDurationABC.output){
  
  cures <- getDurationABC.output$Treated
  years <- getDurationABC.output$Years
  duration.est <- getDurationABC.output$duration
  sim.Prevs <- getDurationABC.output$sim.Prevalences
  sim.Durations <- getDurationABC.output$sim.Durations
  
  sim.Prevs.1 <- sim.Prevs[,1]
  M <- length(sim.Prevs.1)
  
  # Approach 1: vary P_1 and use duration estimate to get prevalence over time
  ann.prev.1 <- matrix(0,nr=M,nc=length(years))
  for(i in 1:M){
    ann.prev.1[i,] <- calcAnnPrev(duration.est,prev1=sim.Prevs.1[i],cures=cures,years=years,unDxTBdur=getDurationABC.output$Duration.u)
  }
  
  # Approach 2: sample over both P_1 and duration
  ann.prev.2 <- matrix(0,nr=M,nc=length(years))
  # randomly select M prevalences and M durations to simulate from
  inds.prev <- sample(c(1:M),replace=T)
  inds.duration <- sample(c(1:M),replace=T)
  
  for(i in 1:M){
    ann.prev.2[i,] <- calcAnnPrev(duration=sim.Durations[inds.duration[i]],prev1=sim.Prevs.1[inds.prev[i]],
                                  cures=cures,years=years,unDxTBdur=getDurationABC.output$Duration.u)
  }
  
  return(list(ann.prev.1=ann.prev.1,ann.prev.2=ann.prev.2))
  
}

# function to convert prevalence data from CalcAnnPrev to long format
getPrevLong <- function(prev.dat.wide,y1,y2,pop.dat,country.name){
  
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
  
  return(long.prev.dat.pop)
}

getPrevCIs <- function(prev.ests,prev.dat,cures,y1,y2,pop.dat,dur.ests,unDxTBdur=5){
  # prev.dat is the data from getPrevLong
  # prev.ests is data from getSummaryData
  # dur.ests are the estimated durations from the main simulation
  
  prev.ests.2 <- prev.ests %>% 
    filter(`Assumed duration of untreated disease`==unDxTBdur) %>%
    select(-c(`Assumed duration of untreated disease`,`Duration of treated disease`))
  
  prev.CIs <- data.frame(t(apply(prev.ests.2,2,quantile,c(0.025,0.975))))
  names(prev.CIs) <- c("LCL","UCL")
  prev.mean <- filter(prev.dat,`Untreated duration`==unDxTBdur)
  
  cures.frame <- data.frame(c(cures,NA))
  names(cures.frame) <- c("Treated")
  
  prev.all <- bind_cols(prev.mean,prev.CIs)
  prev.all <- bind_cols(prev.all,cures.frame)
  
  prev.all$`Proportion treated` <- prev.all$Treated/prev.all$Prevalence
  
  prev.all$Prev100k <- 10^8*(prev.all$Prevalence/pop.dat$Total.Population[pop.dat$label %in% c(y1:y2)])
  prev.all$UCL.100k <- 10^8*(prev.all$UCL/pop.dat$Total.Population[pop.dat$label %in% c(y1:y2)])
  prev.all$LCL.100k <- 10^8*(prev.all$LCL/pop.dat$Total.Population[pop.dat$label %in% c(y1:y2)])
  prev.all$Incidence100k <- prev.all$Prev100k/dur.ests$`Duration of treated disease`[dur.ests$`Assumed duration of untreated disease`==unDxTBdur] 
  prev.all$PrevDelta <- c(NA,diff(prev.all$Prev100k))
  prev.all$PropIncTreated <- (10^8)*prev.all$Treated/(prev.all$Incidence100k*pop.dat$Total.Population[pop.dat$label %in% c(y1:y2)])
  
  
  return(prev.all)
  
}
