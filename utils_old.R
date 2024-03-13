# function to account for uncertainty in first and second prevalence survey values
## uses calcAnnPrev() and samples first prevalence from distn to capture uncertainty in that value
## then it assigns a weight to each potential duration based on the pdf of the second prevalence value
getDurationPDF <- function(country,pop.data,duration.range=seq(1,6,0.01),unDxTBdur=5,M=1000){
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
  prev.1.tmps <- rnorm(M,mean=dat.prev.tmp$`1st Prevalence`,sd=sd.1) # store the simulated prevalence values (per 100k)
  cures.tmp <- dat.cure$cured[dat.cure$Country==country] # number of cures for the country (in thousands)
  years.tmp <- dat.cure$Year[dat.cure$Country==country] # years of cure data
  weights.d <-  r.weights <- matrix(0,nr=M,nc=length(duration.range)) #space to store weights for picking durations
  dur.est.sim <- NULL #vector to store duration selected for each simulation
  for(i in 1:M){
    prev.1 <- prev.1.tmps[i]
    prev.1.n <- (prev.1/100000)*(pop1/1000) #total number of prevalent cases in y1
    
    # get list with each item being a vector of prevalences over the years between surveys
    prevs.tmp.n <- lapply(duration.range,calcAnnPrev,prev1=prev.1.n,cures=cures.tmp,years=years.tmp,unDxTBdur=unDxTBdur)
    prevs.tmp.n.2 <- matrix(unlist(prevs.tmp.n),nr=length(duration.range),byrow=T)
    
    # now I just need the last estimated prevalence for each of these
    prev.2.n <- prevs.tmp.n.2[,dim(prevs.tmp.n.2)[2]]
    
    # need to get these prevalences back to being population adjusted to match distn above
    prev.2.tmp <- (1000*prev.2.n/pop2)*100000
    
    # get probability of each of these being observed, given normal distribution of P_2
    weights.d[i,] <- dnorm(prev.2.tmp,mean=dat.prev.tmp$`2nd Survey Prevalence`,sd=sd.2)
    
    # reweight to sum to 1
    r.weights[i,] <- weights.d[i,]/sum(weights.d[i,])
    
    # get best duration estimate for this simulation
    dur.est.sim[i] <- duration.range[which.max(r.weights[i,])]
    ###note: probably want to change this to accept/reject durations, 
    ####rather than pick best from each sim
  }
  
  # now get estimate of duration and empirical CIs
  dur.est <- mean(dur.est.sim)
  CIs <- quantile(dur.est.sim,c(0.025,0.975))
  
  # get prevalences for each of these duration estimates (estimate and CI values)
  prevs.est <- calcAnnPrev(duration=dur.est,prev1=prev1.n,cures=cures.tmp,
                           years=years.tmp,unDxTBdur=unDxTBdur)
  
  ### fix this!!! Makes all start at observed prevalence ###
  prevs.ucl <- calcAnnPrev(duration=CIs[1],prev1=prev1.n,cures=cures.tmp,
                           years=years.tmp,unDxTBdur=unDxTBdur)
  prevs.lcl <- calcAnnPrev(duration=CIs[2],prev1=prev1.n,cures=cures.tmp,
                           years=years.tmp,unDxTBdur=unDxTBdur)
  ### these are number of cases/1000
  
  # estimated incidence/1000
  ann.inc <- prevs.est*(1/dur.est)
  
  # add in prev and inc per 100k population
  pops <- pop.data$Total.Population[which(pop.data$label==y1):which(pop.data$label==y2)]
  prev.100k  <- (prevs.est/pops)*(10^8)
  prev.100k.lcl  <- (prevs.lcl/pops)*(10^8)
  prev.100k.ucl  <- (prevs.ucl/pops)*(10^8)
  inc.100k <- (ann.inc/pops)*(10^8)
  
  ### NOTE: return all the prevs that were estimated in each simulation ###  
  # data frame with all results
  l <- length(prevs.est)
  results.frame <- data.frame(Country=rep(country,l),Year=seq(y1,y2,1),Prevalence=prevs.est,Prevalence.lcl=prevs.lcl,Prevalence.ucl=prevs.ucl,
                              Treated=c(cures.tmp,NA),Incidence=ann.inc,Duration=rep(dur.est,l),
                              Prev.pop=prev.100k,Prev.pop.lcl=prev.100k.lcl,Prev.pop.ucl=prev.100k.ucl,Inc.pop=inc.100k,
                              Duration.u=rep(unDxTBdur,l),Cure.Prev=c(cures.tmp/prevs.est[1:(l-1)],NA),Cure.Inc=c(cures.tmp/ann.inc[1:(l-1)],NA))
  
  return(list(all.results=results.frame,dur.ests=dur.est,CIs=CIs,sim.durations=dur.est.sim,sim.weights=r.weights,sim.prev1=prev.1.tmps))
}
#######################################
# function to get the duration, incidences, annual prevalences, etc for a country
getDuration2 <- function(country.tmp,pop.data,unDxTBdur=5){
  
  # get years of surveys
  y1 <- dat.prev$`First prevalence survey year`[dat.prev$Country==country.tmp]
  y2 <- dat.prev$`Second prevalence survey year`[dat.prev$Country==country.tmp]
  
  # get first and second prevalence survey estimates (per 100k)
  prev1 <- dat.prev$`1st Prevalence`[dat.prev$Country==country.tmp]
  prev2 <- dat.prev$`2nd Survey Prevalence`[dat.prev$Country==country.tmp]
  
  # convert to number of cases in thousands
  prev1.n <- (prev1/100000)*pop.data$Total.Population[pop.data$label==y1]/1000
  prev2.n <- (prev2/100000)*pop.data$Total.Population[pop.data$label==y2]/1000
  
  # find the duration of disease (i.e. ratio of prevalence to incidence) that gets us closest to the second prevalence numbers
  duration.test <- seq(1,5,0.01)
  
  treated <- dat.cure$cured[dat.cure$Country==country.tmp]
  years.tmp <- dat.cure$Year[dat.cure$Country==country.tmp]
  
  dist.tmp <- prev2.tmp <- NULL
  for(i in 1:length(duration.test)){
    tmp.prevs <- calcAnnPrev(duration.test[i],prev1.n,treated,years=years.tmp,unDxTBdur=unDxTBdur)
    l <- length(tmp.prevs)
    prev2.tmp[i] <- tmp.prevs[l]
    dist.tmp[i] <- abs(prev2.tmp[i]-prev2.n)
  }
  
  dur.est <- duration.test[which.min(dist.tmp)]
  ann.prevs <- calcAnnPrev(dur.est,prev1.n,dat.cure$cured[dat.cure$Country==country.tmp],years=dat.cure$Year[dat.cure$Country==country.tmp],unDxTBdur=unDxTBdur)
  ann.inc <- ann.prevs*(1/dur.est)
  
  # add in prev and inc per 100k population
  pops <- pop.data$Total.Population[which(pop.data$label==y1):which(pop.data$label==y2)]
  prev.100k  <- (ann.prevs/pops)*(10^8)
  inc.100k <- (ann.inc/pops)*(10^8)
  
  results.frame <- data.frame(Country=rep(country.tmp,l),Year=seq(y1,y2,1),Prevalence=ann.prevs,Treated=c(treated,NA),Incidence=ann.inc,Duration=rep(dur.est,l),
                              Prev.pop=prev.100k,Inc.pop=inc.100k,
                              Duration.u=rep(unDxTBdur,l),Cure.Prev=c(treated/ann.prevs[1:(l-1)],NA),Cure.Inc=c(treated/ann.inc[1:(l-1)],NA))
  
  return(results.frame)
  
}

################################################

getDuration <- function(r,prev2,tot.pop,cure.rate,num.deaths,notifications,nat.recovery=0){
  # adds self cure-nat.recovery is taken from entire prevalent pool of cases
  # takes the slope as an input
  # note that prev are rates per 100k
  # time1, time2 are the dates of the prevalence surveys
  # tot.pop: total population of the country
  # cure.rate: proportion of notified TB cases that complete treatment
  # num.deaths: raw number of deaths attributed to TB
  # notifications: annual number of notified cases of TB
  
  # now get the estimated prevalence the year before the last prev survey
  prev.new <- prev2/r
  
  # convert this to the total number of cases
  N.2 <- prev2*tot.pop/100000
  N.new <- prev.new*tot.pop/100000
  
  # number of treatment completions
  cured <- notifications*cure.rate
  
  # number of natural cures from the prevalent cases
  if(nat.recovery > 0){ 
    nat.cures <- nat.recovery*N.new
  } else nat.cures <- 0
  
  # estimated incidence
  Inc.cases <- N.2-(N.new-cured-num.deaths-nat.cures)
  
  # incidence rate
  Inc <- 100000*Inc.cases/tot.pop
  
  # duration
  duration <- prev.new/Inc
  
  # removal rate (divided by prevalence)
  removal.rate.p <- (num.deaths+cured+nat.cures)/N.2
  
  # removal rate (divided by incidence)
  removal.rate.i <- (num.deaths+cured+nat.cures)/Inc.cases
  
  return(list(rate=r,prev2=prev2,Num.cases.2=N.2,Num.cases.prior=N.new,N.cured=cured,incident.cases=Inc.cases,
              incident.rate=Inc,duration=duration,prevalent.cases.target=N.new,
              removal.rate.p=removal.rate.p,removal.rate.i=removal.rate.i,cure.rate=cure.rate,num.deaths=num.deaths,notifications=notifications,
              tot.pop=tot.pop))
}
# Function to calculate CIs based on uncertainty in the prevalence estimates
calcCIs <- function(r,prev2,tot.pop,cure.rate,num.deaths,notifications,prev.lcl,prev.ucl,nat.recovery=0){
  # this function only considers uncertainty in the second prevalence estimate
  
  # 95% CI seems to be based on log normal distribution of prevalence
  log.mean <- log(prev2)
  log.se <- (log(prev.ucl)-log(prev2))/1.96
  
  tmp.prev <- exp(rnorm(10000,log.mean,log.se))
  
  tmp.results <- list()
  tmp.durations <- tmp.incidences <- NULL
  for(i in 1:length(tmp.prev)){
    tmp.results[[i]] <- getDuration(r,tmp.prev[i],tot.pop,cure.rate,num.deaths,notifications,nat.recovery=nat.recovery)
    tmp.durations[i] <- tmp.results[[i]]$duration
    tmp.incidences[i] <- tmp.results[[i]]$incident.rate
  }
  duration.lcl <- quantile(tmp.durations,0.025)
  duration.ucl <- quantile(tmp.durations,0.975)
  
  incidence.lcl <- quantile(tmp.incidences,0.025)
  incidence.ucl <- quantile(tmp.incidences,0.975)
  
  return(list(duration.lcl=duration.lcl,duration.ucl=duration.ucl,
              incidence.lcl=incidence.lcl,incidence.ucl=incidence.ucl,all.sims=tmp.results))
}

# function to combine results into a vector
combineResults <- function(getDurationResults,calcCIsResults,countryName,prev1){
  #function takes results from main functions to calculate results
  
  res.vec <- round(c(getDurationResults$tot.pop,prev1,getDurationResults$prev2,getDurationResults$Num.cases.2,
                     getDurationResults$Num.cases.prior,getDurationResults$rate,getDurationResults$notifications,
                     getDurationResults$cure.rate,
                     getDurationResults$N.cured,getDurationResults$num.deaths,getDurationResults$incident.cases,
                     getDurationResults$incident.rate,calcCIsResults$incidence.lcl,
                     calcCIsResults$incidence.ucl,getDurationResults$duration,calcCIsResults$duration.lcl,
                     calcCIsResults$duration.ucl,
                     getDurationResults$removal.rate.p,getDurationResults$removal.rate.i),digits=2)
  
  res.vec2 <- c(countryName,res.vec)
  return(res.vec2)
  
}

# older function that calculates the linear slope
getDuration2.old <- function(prev1,prev2,time1,time2,tot.pop,cure.rate,num.deaths,notifications){
  # calculates the slope
  # note that prev are rates per 100k
  # time1, time2 are the dates of the prevalence surveys
  # tot.pop: total population of the country
  # cure.rate: proportion of notified TB cases that complete treatment
  # num.deaths: raw number of deaths attributed to TB
  # notifications: annual number of notified cases of TB
  
  # first get slope describing change in prevalence over time
  m <- (prev2-prev1)/(time2-time1)
  
  # now get the estimated prevalence the year before the last prev survey
  prev.new <- prev2-m
  
  # convert this to the total number of cases
  N.2 <- prev2*tot.pop/100000
  N.new <- prev.new*tot.pop/100000
  
  # number of treatment completions
  cured <- notifications*cure.rate
  
  # estimated incidence
  Inc.cases <- N.2-(N.new-cured-num.deaths)
  
  # incidence rate
  Inc <- 100000*Inc.cases/tot.pop
  
  # duration
  duration <- prev.new/Inc
  
  return(list(slope=m,incident.cases=Inc.cases,incident.rate=Inc,duration=duration,prevalent.cases.target=N.new))
}

# function to calculate duration
getDuration.old <- function(r,prev2,tot.pop,cure.rate,num.deaths,notifications,nat.recovery=0){
  # takes the slope as an input
  # note that prev are rates per 100k
  # time1, time2 are the dates of the prevalence surveys
  # tot.pop: total population of the country
  # cure.rate: proportion of notified TB cases that complete treatment
  # num.deaths: raw number of deaths attributed to TB
  # notifications: annual number of notified cases of TB
  
  if(nat.recovery!=0){
    cure.rate <- cure.rate+nat.recovery
  }
  
  # now get the estimated prevalence the year before the last prev survey
  prev.new <- prev2/r
  
  # convert this to the total number of cases
  N.2 <- prev2*tot.pop/100000
  N.new <- prev.new*tot.pop/100000
  
  # number of treatment completions
  cured <- notifications*cure.rate
  
  # estimated incidence
  Inc.cases <- N.2-(N.new-cured-num.deaths)
  
  # incidence rate
  Inc <- 100000*Inc.cases/tot.pop
  
  # duration
  duration <- prev.new/Inc
  
  # removal rate (divided by prevalence)
  removal.rate.p <- (num.deaths+cured)/N.2
  
  # removal rate (divided by incidence)
  removal.rate.i <- (num.deaths+cured)/Inc.cases
  
  return(list(rate=r,prev2=prev2,Num.cases.2=N.2,Num.cases.prior=N.new,N.cured=cured,incident.cases=Inc.cases,
              incident.rate=Inc,duration=duration,prevalent.cases.target=N.new,
              removal.rate.p=removal.rate.p,removal.rate.i=removal.rate.i,cure.rate=cure.rate,num.deaths=num.deaths,notifications=notifications,
              tot.pop=tot.pop))
}

