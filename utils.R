
#function to get rate of change assuming multiplicative rate of change
getRateofChange <- function(prev1,prev2,year1,year2){
  time.period <- year2-year1
  r <- (prev2/prev1)^(1/time.period)
  
  return(r)
}

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
getDuration2 <- function(prev1,prev2,time1,time2,tot.pop,cure.rate,num.deaths,notifications){
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

