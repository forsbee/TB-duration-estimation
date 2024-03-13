test <- getDurationPDF(country="Vietnam",dat.prev=dat.prev,M=10000)

wts <- test[[1]]
prevs <- test[[2]]

weights.k <- apply(wts,2,mean)
summary(weights.k)

weights.k <- weights.k/sum(weights.k)
summary(weights.k)
sum(weights.k)
duration.range[which.max(weights.k)]

# now try reweighting rows first
r.weights.k <- matrix(0,nr=10000,nc=dim(wts)[2])
for(i in 1:10000){
  r.weights.k[i,] <- wts[i,]/sum(wts[i,])
}

weights.k.2 <- apply(r.weights.k,2,mean)
summary(weights.k.2)

duration.range[which.max(weights.k.2)]

plot(duration.range,weights.k.2)
points(duration.range,weights.k,col="blue")


# try to determine the 95% CI from this
lower.CI <- duration.range[which.min(abs(cumsum(weights.k)-0.025))]
upper.CI <- duration.range[which.min(abs(cumsum(weights.k)-0.975))]
est <- duration.range[which.min(abs(cumsum(weights.k)-0.5))]
duration.range[which.max(weights.k)]
c(est,lower.CI,upper.CI)

###########################################
# another approach is to get est and then average those (cannot do that so easily with CIs-would need to get quantiles from distribution of ests)
dim(r.weights.k)
inds <- apply(r.weights.k,1,which.max)
length(inds)
ests.sim <- duration.range[inds]
summary(ests.sim)
CIs <- quantile(ests.sim,c(0.025,0.975))
############################################

# alternatively, should I weight the rows of the weight matrix by the probability of the prevalence? probably not since sampling randomly is already undersampling 
## extreme values and this will only further undervalue those-empirical approach above is probably best.
p.prev <- dnorm(prevs,mean=dat.prev.tmp$`1st Prevalence`,sd=sd.1)
length(p.prev)
sum(p.prev)
p.prev.2 <- p.prev/sum(p.prev)
sum(p.prev.2)

weights.k.3 <- t(matrix(p.prev.2,nc=1))%*%r.weights.k
points(duration.range,weights.k.3,col="red")
duration.range[which.max(weights.k.3)]

# try to determine the 95% CI from this
lower.CI <- duration.range[which.min(abs(cumsum(weights.k.3)-0.025))]
upper.CI <- duration.range[which.min(abs(cumsum(weights.k.3)-0.975))]
est <- duration.range[which.min(abs(cumsum(weights.k.3)-0.5))]

c(est,lower.CI,upper.CI)
