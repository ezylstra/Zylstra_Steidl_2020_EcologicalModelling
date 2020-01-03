######################################################################################
#Zylstra, ER, and RJ Steidl. 2020. A Bayesian state-space model for seasonal growth 
#of terrestrial vertebrates. Ecological Modelling.

#Simulating seasonal growth and estimating parameters from seasonal and 
#non-seasonal von Bertalanffy models        
######################################################################################

#------------------------------------------------------------------------------------------# 
# Note that notation below is different than that used in the manuscript.
# Here, asymptotic length is represented by Linf (in the manuscript, we used L)
# Here, length at first capture is represented by L1 (in the manscript, we use lambda)
#------------------------------------------------------------------------------------------#


# install.packages(c('plyr','''reshape2','coda','jagsUI'))
library(plyr)
library(reshape2)
library(coda)
library(jagsUI)

#---------------- Survey design parameters -----------------------------------------#
#------ 16 total occasions, spanning 2-year period (2014-2015)
#------ Captures in 8 consective months per year (Mar - Oct)

meanints.m <- c(rep(1,7),5,rep(1,7))
meanints.y <- meanints.m/12
occdates <- as.Date(c('2014-03-01','2014-04-01','2014-05-01','2014-06-01',
                      '2014-07-01','2014-08-01','2014-09-01','2014-10-01',
                      '2015-03-01','2015-04-01','2015-05-01','2015-06-01',
                      '2015-07-01','2015-08-01','2015-09-01','2015-10-01'))
#Fraction of year associated with each capture occasion
tyear <- as.numeric(format(occdates,'%j'))/365

#---------------- Growth parameters -----------------------------------------#
#------ Individual heterogeneity in asymptotic length but not k

Linf.mu <- 50        #mean asymptotic size
Linf.sd <- 2         #individual heterogeneity in asymptotic size
L1.mu <- 40          #mean length at first capture
L1.sd <- 4           #individual heterogeneity in first capture
k <- 3               #growth parameter
C <- 1               #1 = seasonal cessation of growth; 0 = constant growth (no seasonal variation)
ts <- 0.60           #time of peak growth (in fractions of a year)
meas.sd <- 1.2       #measurement error

#---------------- Other simulation parameters -------------------------------#

#Number of individuals
  ncap.occ <- 40  #Number of frogs first captured during each occasion
  nocc <- length(occdates)
  nind <- (nocc-1)*ncap.occ
  firstocc <- rep(1:(nocc-1),each=ncap.occ)

#Survival probabilities
  #50% annual survival
  #Assume individuals survive first interval with probability=1 (need at least one recap)
  #Assume survival during all other intervals (winter or not) the same
  annprob <- 0.5
  monthlyprob <- annprob^(1/7)   #there are 8 intervals per year, so 7 with survival <1
  monthlyprobsurvival <- c(1,rep(monthlyprob,nocc-2))

#Recapture probabilities
  #Constant 70% recapture probability [after initial recapture]
  probnextcap <- c(1,rep(0.7,nocc-2))  

#---------------- JAGS model for seasonal growth, estimating C
sink("VBgrowth_Season_sims.txt")
cat("
  model{

  for(m in 1:nmeas){
    y[m] ~ dnorm(mu[m],meas.tau)
    mu[m] <- L1[ind[m]] + (Linf[ind[m]]-L1[ind[m]])*(1-exp(-k*d[m] -
             C*k/(2*pi)*sin(2*pi*(t[m]+d[m]-ts)) + C*k/(2*pi)*sin(2*pi*(t[m]-ts))))

  }

  for(i in 1:nind){
    L1[i] ~ dnorm(L1.mu,L1.tau)T(0,Linf[i])
    Linf[i] ~ dnorm(Linf.mu,Linf.tau)T(0,)
  }

  log(k) <- logk
  L1.mu ~ dnorm(0,0.001)
  Linf.mu ~ dnorm(0,0.001)
  logk ~ dnorm(0,0.001)

  meas.tau <- 1/(meas.sd*meas.sd)
  meas.sd ~ dt(0,pow(5,-2),1)T(0,)

  L1.tau <- 1/(L1.sd*L1.sd)
  L1.sd ~ dt(0,pow(5,-2),1)T(0,)

  Linf.tau <- 1/(Linf.sd*Linf.sd)
  Linf.sd ~ dt(0,pow(5,-2),1)T(0,)

  pi <- 3.14159265359
  ts ~ dunif(0,1)
  C ~ dunif(0,1)

  #derived parameters:
  ts.julian <- ts*365

  } #model
  ",fill=TRUE)
sink()

#---------------- JAGS model for seasonal growth, fixed C
sink("VBgrowth_SeasonFixedC_sims.txt")
cat("
  model{

  for(m in 1:nmeas){
    y[m] ~ dnorm(mu[m],meas.tau)
    mu[m] <- L1[ind[m]] + (Linf[ind[m]]-L1[ind[m]])*(1-exp(-k*d[m] -
             C*k/(2*pi)*sin(2*pi*(t[m]+d[m]-ts)) + C*k/(2*pi)*sin(2*pi*(t[m]-ts))))

  }

  for(i in 1:nind){
    L1[i] ~ dnorm(L1.mu,L1.tau)T(0,Linf[i])
    Linf[i] ~ dnorm(Linf.mu,Linf.tau)T(0,)
  }

  log(k) <- logk
  L1.mu ~ dnorm(0,0.001)
  Linf.mu ~ dnorm(0,0.001)
  logk ~ dnorm(0,0.001)

  meas.tau <- 1/(meas.sd*meas.sd)
  meas.sd ~ dt(0,pow(5,-2),1)T(0,)

  L1.tau <- 1/(L1.sd*L1.sd)
  L1.sd ~ dt(0,pow(5,-2),1)T(0,)

  Linf.tau <- 1/(Linf.sd*Linf.sd)
  Linf.sd ~ dt(0,pow(5,-2),1)T(0,)

  pi <- 3.14159265359
  ts ~ dunif(0,1)
  C <- 1

  #derived parameters:
  ts.julian <- ts*365

  } #model
  ",fill=TRUE)
sink()
  
#---------------- JAGS model for seasonal growth, fixed C and no individual heterogeneity in L
sink("VBgrowth_SeasonFixedCConstL_sims.txt")
cat("
  model{

  for(m in 1:nmeas){
    y[m] ~ dnorm(mu[m],meas.tau)
    mu[m] <- L1[ind[m]] + (Linf.mu-L1[ind[m]])*(1-exp(-k*d[m] -
             C*k/(2*pi)*sin(2*pi*(t[m]+d[m]-ts)) + C*k/(2*pi)*sin(2*pi*(t[m]-ts))))

  }

  for(i in 1:nind){
    L1[i] ~ dnorm(L1.mu,L1.tau)T(0,Linf.mu)
  }

  log(k) <- logk
  L1.mu ~ dnorm(0,0.001)
  Linf.mu ~ dnorm(0,0.001)
  logk ~ dnorm(0,0.001)

  meas.tau <- 1/(meas.sd*meas.sd)
  meas.sd ~ dt(0,pow(5,-2),1)T(0,)

  L1.tau <- 1/(L1.sd*L1.sd)
  L1.sd ~ dt(0,pow(5,-2),1)T(0,)

  pi <- 3.14159265359
  ts ~ dunif(0,1)
  C <- 1

  #derived parameters:
  ts.julian <- ts*365

  } #model
  ",fill=TRUE)
sink()
 
#---------------- JAGS model for no seasonal variation in growth
sink("VBgrowth_NoSeason_sims.txt")
cat("
  model{

  for(m in 1:nmeas){
    y[m] ~ dnorm(mu[m],meas.tau)
    mu[m] <- L1[ind[m]] + (Linf[ind[m]]-L1[ind[m]])*(1-exp(-k*d[m]))

  }

  for(i in 1:nind){
    L1[i] ~ dnorm(L1.mu,L1.tau)T(0,Linf[i])
    Linf[i] ~ dnorm(Linf.mu,Linf.tau)T(0,)
  }

  log(k) <- logk
  L1.mu ~ dnorm(0,0.001)
  Linf.mu ~ dnorm(0,0.001)
  logk ~ dnorm(0,0.001)

  meas.tau <- 1/(meas.sd*meas.sd)
  meas.sd ~ dt(0,pow(2.5,-2),1)T(0,)

  L1.tau <- 1/(L1.sd*L1.sd)
  L1.sd ~ dt(0,pow(10,-2),1)T(0,)

  Linf.tau <- 1/(Linf.sd*Linf.sd)
  Linf.sd ~ dt(0,pow(2.5,-2),1)T(0,)

  } #model
  ",fill=TRUE)
sink()

#---------------- parameters/settings for JAGS ------------------------------#
inits.seasEstC <- function(){list(Linf.mu = runif(1,45,55),
                                  Linf.sd = runif(1,0,3),
                                  L1.mu = runif(1,35,48),
                                  L1.sd = runif(1,0,4),
                                  logk = runif(1,0,2), 
                                  C = runif(1,0.5,1),
                                  ts = runif(1,0.4,0.8),
                                  meas.sd = runif(1,0,3))}
inits.seasFixC <- function(){list(Linf.mu = runif(1,45,55),
                                  Linf.sd = runif(1,0,3),
                                  L1.mu = runif(1,35,48),
                                  L1.sd = runif(1,0,4),
                                  logk = runif(1,0,2), 
                                  ts = runif(1,0.4,0.8),
                                  meas.sd = runif(1,0,3))}
inits.seasConL <- function(){list(Linf.mu = runif(1,45,55),
                                  L1.mu = runif(1,35,48),
                                  L1.sd = runif(1,0,4),
                                  logk = runif(1,0,2),  
                                  ts = runif(1,0.4,0.8),
                                  meas.sd = runif(1,0,3))}
inits.const <- function(){list(Linf.mu = runif(1,45,55),
                               Linf.sd = runif(1,0,3),
                               L1.mu = runif(1,35,48),
                               L1.sd = runif(1,0,4),
                               logk = runif(1,0,2),  
                               meas.sd = runif(1,0,3))}

#---------------- objects for storing simulation results --------------------#
simresults.seasEstC <- matrix(NA,nrow=1,ncol=31)
simresults.seasEstC <- as.data.frame(simresults.seasEstC)
names(simresults.seasEstC) <- c('sim.no','ni','nt','nb','na',
                                'sim.Linf','sim.Linf.sd',
                                'sim.L1','sim.L1.sd',
                                'sim.k','sim.meas.sd',
                                'sim.C','sim.ts',
                                'Linf.mu','Linf.sd',
                                'L1.mu','L1.sd',
                                'k','meas.sd',
                                'C','ts',
                                'Linf.mu.cov','Linf.sd.cov',
                                'L1.mu.cov','L1.sd.cov',
                                'k.cov','meas.sd.cov',
                                'C.cov','ts.cov',
                                'Rhat.mean','Rhat.max')
params.seasEstC <- c('Linf.mu','Linf.sd','L1.mu','L1.sd','k','meas.sd','C','ts')

simresults.seasFixC <- matrix(NA,nrow=1,ncol=28)
simresults.seasFixC <- as.data.frame(simresults.seasFixC)
names(simresults.seasFixC) <- c('sim.no','ni','nt','nb','na',
                                'sim.Linf','sim.Linf.sd',
                                'sim.L1','sim.L1.sd',
                                'sim.k','sim.meas.sd',
                                'sim.ts',
                                'Linf.mu','Linf.sd',
                                'L1.mu','L1.sd',
                                'k','meas.sd',
                                'ts',
                                'Linf.mu.cov','Linf.sd.cov',
                                'L1.mu.cov','L1.sd.cov',
                                'k.cov','meas.sd.cov',
                                'ts.cov',
                                'Rhat.mean','Rhat.max')
params.seasFixC <- c('Linf.mu','Linf.sd','L1.mu','L1.sd','k','meas.sd','ts')

simresults.seasConL <- matrix(NA,nrow=1,ncol=26)
simresults.seasConL <- as.data.frame(simresults.seasConL)
names(simresults.seasConL) <- c('sim.no','ni','nt','nb','na',
                                'sim.Linf','sim.Linf.sd',
                                'sim.L1','sim.L1.sd',
                                'sim.k','sim.meas.sd',
                                'sim.ts',
                                'Linf.mu',
                                'L1.mu','L1.sd',
                                'k','meas.sd',
                                'ts',
                                'Linf.mu.cov',
                                'L1.mu.cov','L1.sd.cov',
                                'k.cov','meas.sd.cov',
                                'ts.cov',
                                'Rhat.mean','Rhat.max')
params.seasConL <- c('Linf.mu','L1.mu','L1.sd','k','meas.sd','ts')

simresults.const <- matrix(NA,nrow=1,ncol=25)
simresults.const <- as.data.frame(simresults.const)
names(simresults.const) <- c('sim.no','ni','nt','nb','na',
                            'sim.Linf','sim.Linf.sd',
                            'sim.L1','sim.L1.sd',
                            'sim.k','sim.meas.sd',
                            'Linf.mu','Linf.sd',
                            'L1.mu','L1.sd',
                            'k','meas.sd',
                            'Linf.mu.cov','Linf.sd.cov',
                            'L1.mu.cov','L1.sd.cov',
                            'k.cov','meas.sd.cov',
                            'Rhat.mean','Rhat.max')
params.const <- c('Linf.mu','Linf.sd','L1.mu','L1.sd','k','meas.sd')

#---------------- Simulating seasonal growth data --------------------#

no.sims <- 225

#MCMC settings
ni <- 5000; nt <- 5; nc <- 3;  nb <- na <- 1000
ni.tot <- nb + ni

for(s in 1:no.sims){

  siml <- matrix(NA,nrow=nind,ncol=nocc)
  l1 <- round(rnorm(nind,mean=L1.mu,sd=L1.sd))
  t1 <- firstocc
  for(i in 1:nrow(siml)){
    siml[i,t1[i]] <- l1[i]
  }
  sims <- as.data.frame(siml)
  names(sims) <- paste('l',1:ncol(sims),sep='')  
  #Sims is a dataframe with initial SVL for each individual

  #Simulate individual asymptotic size
  simLinf <- rnorm(nind,mean=Linf.mu,sd=Linf.sd)

  #Simulate true lengths of each individual at each capture occasion:
  for(i in 1:nind){
    for(j in (t1[i]+1):nocc){
      sims[i,j] <- l1[i] + (simLinf[i]-l1[i])*(1-exp(-k*(sum(meanints.y[t1[i]:(j-1)])) - 
                   C*k/(2*pi)*sin(2*pi*(tyear[t1[i]]+sum(meanints.y[t1[i]:(j-1)])-ts)) +
                   C*k/(2*pi)*sin(2*pi*(tyear[t1[i]]-ts)))) 
    }
  }
  #Add in measurement error:
  for(i in 1:nind){
    for(j in t1[i]:nocc){
      sims[i,j] <- round(rnorm(1,sims[i,j],meas.sd))
    }
  }
  #Only retain lengths when frogs are alive (and make sure frogs stay dead after alive=0):
  for(i in 1:nind){
    for(j in t1[i]:(nocc-1)){
      alive <- rbinom(1,1,monthlyprobsurvival[j-t1[i]+1])
      sims[i,j+1] <- ifelse(alive==0 | is.na(sims[i,j]),NA,sims[i,j+1])
    }
  }  
  #Only retain lengths when frogs captured:
  for(i in 1:nind){
    for(j in t1[i]:(nocc-1)){
      cap <- rbinom(1,1,probnextcap[j-t1[i]+1])
      sims[i,j+1] <- ifelse(cap==0,NA,sims[i,j+1])
    }
  }  

  sims$ncap <- rowSums(!is.na(sims))
  sims$id <- 1:nrow(sims)

  #Create dataframe with information about the first capture of each frog:
  df.frog <- data.frame(id=sims$id,t1=t1,t=tyear[t1],ncap=sims$ncap)
  
  #Create dataframe with each capture (including first) of each frog:
  df.cap <- melt(sims,id.vars=c('id','ncap'),value.name='svl',variable.name='occ')
  df.cap <- df.cap[with(df.cap,order(id,occ)),]
  df.cap <- df.cap[!is.na(df.cap$svl),]
  df.cap <- join(df.cap,df.frog[,c('id','t1','t')],by='id',type='left')
  df.cap$occ <- as.numeric(substr(df.cap$occ,2,nchar(as.character(df.cap$occ))))
  for(i in 1:nrow(df.cap)){
    df.cap$d[i] <- ifelse(df.cap$occ[i]==df.cap$t1[i],0,sum(meanints.y[df.cap$t1[i]:(df.cap$occ[i]-1)])) 
  }

  nmeas <- nrow(df.cap)
  
  #----- Run seasonal-EstC JAGS model
  seasmodel.data <- list(y=df.cap$svl,d=df.cap$d,t=df.cap$t,ind=df.cap$id,nind=nind,nmeas=nmeas)

  seasm <- jags(data=seasmodel.data, inits=inits.seasEstC, parameters.to.save=params.seasEstC,
                model.file='VBgrowth_Season_sims.txt', n.chains=nc, n.iter=ni.tot,
                n.burnin=nb, n.thin=nt, parallel=T, store.data=T)

    simresults.seasEstC[s,'sim.no'] <- s
    simresults.seasEstC[s,'ni'] <- ni 
    simresults.seasEstC[s,'nt'] <- nt 
    simresults.seasEstC[s,'nb'] <- nb
    simresults.seasEstC[s,'na'] <- na
    simresults.seasEstC[s,6:13] <- c(Linf.mu,Linf.sd,L1.mu,L1.sd,k,meas.sd,C,ts)
    indseas <- match(params.seasEstC,colnames(seasm$samples[[1]]))
    indseas <- indseas[!is.na(indseas)]
    simresults.seasEstC[s,14:21] <- summary(seasm$samples)$statistics[indseas,1]
    for(j in 1:8){
      simresults.seasEstC[s,21+j] <- ifelse(get(params.seasEstC[j])>=summary(seasm$samples)$quantiles[indseas[j],1] & 
                                            get(params.seasEstC[j])<=summary(seasm$samples)$quantiles[indseas[j],5],1,0)
    }
    simresults.seasEstC[s,'Rhat.mean'] <- gelman.diag(seasm$samples)$mpsrf
    simresults.seasEstC[s,'Rhat.max'] <- max(gelman.diag(seasm$samples)$psrf[,'Point est.']) 
    
  #----- Run seasonal-FixC JAGS model
  seasmodelFC.data <- list(y=df.cap$svl,d=df.cap$d,t=df.cap$t,ind=df.cap$id,nind=nind,nmeas=nmeas)

  seasFCm <- jags(data=seasmodelFC.data, inits=inits.seasFixC, parameters.to.save=params.seasFixC,
                  model.file='VBgrowth_SeasonFixedC_sims.txt', n.chains=nc, n.iter=ni.tot,
                  n.burnin=nb, n.thin=nt, parallel=T, store.data=T)

    simresults.seasFixC[s,'sim.no'] <- s
    simresults.seasFixC[s,'ni'] <- ni 
    simresults.seasFixC[s,'nt'] <- nt 
    simresults.seasFixC[s,'nb'] <- nb
    simresults.seasFixC[s,'na'] <- na
    simresults.seasFixC[s,6:12] <- c(Linf.mu,Linf.sd,L1.mu,L1.sd,k,meas.sd,ts)
    indseasFC <- match(params.seasFixC,colnames(seasFCm$samples[[1]]))
    indseasFC <- indseasFC[!is.na(indseasFC)]
    simresults.seasFixC[s,13:19] <- summary(seasFCm$samples)$statistics[indseasFC,1]
    for(j in 1:7){
      simresults.seasFixC[s,19+j] <- ifelse(get(params.seasFixC[j])>=summary(seasFCm$samples)$quantiles[indseasFC[j],1] & 
                                            get(params.seasFixC[j])<=summary(seasFCm$samples)$quantiles[indseasFC[j],5],1,0)
    }
    simresults.seasFixC[s,'Rhat.mean'] <- gelman.diag(seasFCm$samples)$mpsrf
    simresults.seasFixC[s,'Rhat.max'] <- max(gelman.diag(seasFCm$samples)$psrf[,'Point est.'])  

  #----- Run seasonal-ConstantL JAGS model
  seasmodelCL.data <- list(y=df.cap$svl,d=df.cap$d,t=df.cap$t,ind=df.cap$id,nind=nind,nmeas=nmeas)

  seasCLm <- jags(data=seasmodelCL.data, inits=inits.seasConL, parameters.to.save=params.seasConL,
                  model.file='VBgrowth_SeasonFixedCConstL_sims.txt', n.chains=nc, n.iter=ni.tot,
                  n.burnin=nb, n.thin=nt, parallel=T, store.data=T)

    simresults.seasConL[s,'sim.no'] <- s
    simresults.seasConL[s,'ni'] <- ni 
    simresults.seasConL[s,'nt'] <- nt 
    simresults.seasConL[s,'nb'] <- nb
    simresults.seasConL[s,'na'] <- na
    simresults.seasConL[s,6:12] <- c(Linf.mu,Linf.sd,L1.mu,L1.sd,k,meas.sd,ts)
    indseasCL <- match(params.seasConL,colnames(seasCLm$samples[[1]]))
    indseasCL <- indseasCL[!is.na(indseasCL)]
    simresults.seasConL[s,13:18] <- summary(seasCLm$samples)$statistics[indseasCL,1]
    for(j in 1:6){
      simresults.seasConL[s,18+j] <- ifelse(get(params.seasConL[j])>=summary(seasCLm$samples)$quantiles[indseasCL[j],1] & 
                                            get(params.seasConL[j])<=summary(seasCLm$samples)$quantiles[indseasCL[j],5],1,0)
    }
    simresults.seasConL[s,'Rhat.mean'] <- gelman.diag(seasCLm$samples)$mpsrf
    simresults.seasConL[s,'Rhat.max'] <- max(gelman.diag(seasCLm$samples)$psrf[,'Point est.'])  
    
  #----- Run non-seasonal JAGS model
  noseasmodel.data <- list(y=df.cap$svl,d=df.cap$d,ind=df.cap$id,nind=nind,nmeas=nmeas)

  noseasm <- jags(data=noseasmodel.data, inits=inits.const, parameters.to.save=params.const,
                model.file='VBgrowth_NoSeason_sims.txt', n.chains=nc, n.iter=ni.tot,
                n.burnin=nb, n.thin=nt, parallel=T, store.data=T)

    simresults.const[s,'sim.no'] <- s
    simresults.const[s,'ni'] <- ni 
    simresults.const[s,'nt'] <- nt 
    simresults.const[s,'nb'] <- nb
    simresults.const[s,'na'] <- na
    simresults.const[s,6:11] <- c(Linf.mu,Linf.sd,L1.mu,L1.sd,k,meas.sd)
    indconst <- match(params.const,colnames(noseasm$samples[[1]]))
    indconst <- indconst[!is.na(indconst)]
    simresults.const[s,12:17] <- summary(noseasm$samples)$statistics[indconst,1]
    for(j in 1:6){
      simresults.const[s,17+j] <- ifelse(get(params.const[j])>=summary(noseasm$samples)$quantiles[indconst[j],1] & 
                                         get(params.const[j])<=summary(noseasm$samples)$quantiles[indconst[j],5],1,0)
    }
    simresults.const[s,'Rhat.mean'] <- gelman.diag(noseasm$samples)$mpsrf
    simresults.const[s,'Rhat.max'] <- max(gelman.diag(noseasm$samples)$psrf[,'Point est.'])  
      
} #s


simresults.const1 <- simresults.const
simresults.seasE1 <- simresults.seasEstC
simresults.seasF1 <- simresults.seasFixC
simresults.seasL1 <- simresults.seasConL

#---------------- Assessing convergence --------------------#
  test1 <- data.frame(seasE=simresults.seasE1$Rhat.mean,const=simresults.const1$Rhat.mean,
                      seasF=simresults.seasF1$Rhat.mean,seasL=simresults.seasL1$Rhat.mean)
  test1$bade <- as.numeric(test1$seasE>1.05)
  test1$badc <- as.numeric(test1$const>1.05)
  test1$badf <- as.numeric(test1$seasF>1.05)
  test1$badl <- as.numeric(test1$seasL>1.05)
  delete <- which(rowSums(test1[,c('bade','badc','badf','badl')])>0)
  
  #Use only those simulated datasets where both models had Rhat<=1.05
  results.c1 <- head(simresults.const1[-delete,],100)
  results.e1 <- head(simresults.seasE1[-delete,],100)
  results.f1 <- head(simresults.seasF1[-delete,],100)
  results.l1 <- head(simresults.seasL1[-delete,],100)
  
#---------------- Evaluating performance of seasonal model, estimating C --------------------#
  #PRB = Percent relative bias = (mean(Est)-True)/True*100
  (prb.e1 <- (colMeans(results.e1[,c('Linf.mu','Linf.sd','L1.mu','L1.sd','k','meas.sd','C','ts')]) -
               results.e1[1,c('sim.Linf','sim.Linf.sd','sim.L1','sim.L1.sd','sim.k','sim.meas.sd','sim.C','sim.ts')])/
               results.e1[1,c('sim.Linf','sim.Linf.sd','sim.L1','sim.L1.sd','sim.k','sim.meas.sd','sim.C','sim.ts')]*100)

  #Coverage probabilities
  colMeans(results.e1[,grep('cov',names(results.e1))])

#---------------- Evaluating performance of seasonal model, fixing C=1 --------------------#
  #PRB = Percent relative bias = (mean(Est)-True)/True*100
  (prb.f1 <- (colMeans(results.f1[,c('Linf.mu','Linf.sd','L1.mu','L1.sd','k','meas.sd','ts')]) -
               results.f1[1,c('sim.Linf','sim.Linf.sd','sim.L1','sim.L1.sd','sim.k','sim.meas.sd','sim.ts')])/
               results.f1[1,c('sim.Linf','sim.Linf.sd','sim.L1','sim.L1.sd','sim.k','sim.meas.sd','sim.ts')]*100)

  #Coverage probabilities
  colMeans(results.f1[,grep('cov',names(results.f1))])

#---------------- Evaluating performance of seasonal model, fixed C=1 and assuming no het in L --------------------#
  #PRB = Percent relative bias = (mean(Est)-True)/True*100
  (prb2.l1 <- (colMeans(results.l1[,c('Linf.mu','L1.mu','L1.sd','k','meas.sd','ts')]) -
               results.l1[1,c('sim.Linf','sim.L1','sim.L1.sd','sim.k','sim.meas.sd','sim.ts')])/
               results.l1[1,c('sim.Linf','sim.L1','sim.L1.sd','sim.k','sim.meas.sd','sim.ts')]*100)

  #Coverage probabilities
  colMeans(results.l1[,grep('cov',names(results.l1))])

#---------------- Evaluating performance of non-seasonal model --------------------#
  #PRB = Percent relative bias = (mean(Est)-True)/True*100
  (prb2.c1 <- (colMeans(results.c1[,c('Linf.mu','Linf.sd','L1.mu','L1.sd','k','meas.sd')]) -
               results.c1[1,c('sim.Linf','sim.Linf.sd','sim.L1','sim.L1.sd','sim.k','sim.meas.sd')])/
               results.c1[1,c('sim.Linf','sim.Linf.sd','sim.L1','sim.L1.sd','sim.k','sim.meas.sd')]*100)

  #Coverage probabilities
  colMeans(results.c1[,grep('cov',names(results.c1))])

#---------------- Histograms ------------------------#
  #Histograms of posterior means for full seasonal model:
  lwidth <- 1.8
  par(mfrow=c(2,3),mar=c(4.5,2,0.5,1),las=1,oma=c(0,3,0,0),cex.lab=1.2,mgp=c(2.8,1,0))
  hist(results.e1$Linf.mu,breaks=seq(49.5,50.5,length=15),col='gray80',main='',xlab=expression(bar(L)),ylab='')
  abline(v=results.e1$sim.Linf,col='blue',lwd=lwidth,lty=2)
  hist(results.e1$k,breaks=seq(2.7,3.3,length=15),col='gray80',main='',xlab='k',ylab='')
  abline(v=results.e1$sim.k,col='blue',lwd=lwidth,lty=2)
  hist(results.e1$C,breaks=seq(0.80,1,length=15),col='gray80',main='',xlab='C',ylab='')
  abline(v=results.e1$sim.C,col='blue',lwd=lwidth,lty=2)
  hist(results.e1$ts,breaks=seq(0.58,0.62,length=15),col='gray80',main='',xlab=expression(paste("t" [s])),ylab='')
  abline(v=results.e1$sim.ts,col='blue',lwd=lwidth,lty=2)
  hist(results.e1$Linf.sd,breaks=seq(1.6,2.4,length=15),col='gray80',main='',xlab=expression(sigma [beta]),ylab='')
  abline(v=results.e1$sim.Linf.sd,col='blue',lwd=lwidth,lty=2)
  hist(results.e1$meas.sd,breaks=seq(1.15,1.35,length=15),col='gray80',main='',xlab=expression(sigma ['y']),ylab='')
  abline(v=results.e1$sim.meas.sd,col='blue',lwd=lwidth,lty=2)
  mtext(text='Frequency',side=2,outer=TRUE,line=1,at=c(0.3,0.8),las=0)

#--------------- Plot predicted growth curves based on model estimates -------------------------#
  Linf <- results.e1$sim.Linf[1]
  k.s <- mean(results.e1$k)
  k.c <- mean(results.c1$k)
  ts <- results.e1$sim.ts[1]
  C <- 1
  int <- 1/365
  Lt1 <- 22.5  #size at metamorphosis
  
  #Indices for day of year
  yrdays <- c(rep(1:365,2),1)  #2013 through 1 Jan 2015
  l13 <- length(yrdays)
  #Set metamorphosis date:
  stdate.s13 <- as.Date('2013-07-01')
  stdate.f13 <- as.Date('2013-08-01')  
  #Set end date:
  enddate <- as.Date('2015-01-01')
  #Convert metamorphosis date to fraction of a year
  day1.s13 <- as.numeric(strftime(stdate.s13,format='%j')); t1.s13 <- day1.s13/365
  day1.f13 <- as.numeric(strftime(stdate.f13,format='%j')); t1.f13 <- day1.f13/365
  #Create vectors for length of individuals in each cohort over time
  L.s13S <- c(Lt1,rep(NA,l13-day1.s13))
  L.f13S <- c(Lt1,rep(NA,l13-day1.f13))
  L.s13C <- c(Lt1,rep(NA,l13-day1.s13))
  L.f13C <- c(Lt1,rep(NA,l13-day1.f13))

  #Predict length for each day after metamorphosis based on estimates from seasonal model:
  for(i in 2:(l13-day1.s13)){  
    L.s13S[i] <- Lt1 + (Linf-Lt1)*(1-exp(-k.s*int*(i-1) - C*k.s/(2*pi)*sin(2*pi*(t1.s13+int*(i-1)-ts)) +
                  C*k.s/(2*pi)*sin(2*pi*(t1.s13-ts)))) 
  }
  for(i in 2:(l13-day1.f13)){  
    L.f13S[i] <- Lt1 + (Linf-Lt1)*(1-exp(-k.s*int*(i-1) - C*k.s/(2*pi)*sin(2*pi*(t1.f13+int*(i-1)-ts)) +
                  C*k.s/(2*pi)*sin(2*pi*(t1.f13-ts)))) 
  }
  #Predict length for each day after metamorphosis based on estimates from non-seasonal model:
  for(i in 2:(l13-day1.s13)){  
    L.s13C[i] <- Lt1 + (Linf-Lt1)*(1-exp(-k.c*int*(i-1)))
  }
  for(i in 2:(l13-day1.f13)){  
    L.f13C[i] <- Lt1 + (Linf-Lt1)*(1-exp(-k.c*int*(i-1)))
  }
  #Dates for plotting
  dateplot.s13 <- seq(from=stdate.s13,to=enddate,by='days')
  dateplot.f13 <- seq(from=stdate.f13,to=enddate,by='days')

#---Plot difference in time to sexual maturity (40mm) and 90% of asymptotic size -- 1 Aug emergence:
  par(mfrow=c(1,1),mar=c(2,2.2,0.5,0.5)+0.1,oma=c(0,0,0,0))
  plot(1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',
       xlim=c(dateplot.f13[1]-10,enddate),ylim=c(21.5,50),bty='n',xaxs='i',yaxs='i')
    usr <- par('usr')  #these are plotting limits (incl extra bit)
    axis(1,at=c(usr[1],usr[2]),tck=F,labels=F)
    axis.Date(1,at=seq(stdate.f13,enddate+1,by='quarter'),labels=seq(2,19,by=3)-1,
              tcl=-0.2,cex.axis=0.75,mgp=c(1.5,0.1,0))
    axis(2,at=c(usr[3],usr[4]),tck=F,labels=F)
    axis(2,at=seq(25,50,by=5),labels=seq(25,50,by=5),las=1,
         tcl=-0.2,cex.axis=0.75,mgp=c(1.5,0.4,0))
    lines(L.f13C~dateplot.f13,lwd=1,col='darkgray')
    lines(L.f13S~dateplot.f13,lwd=1)
    mtext('Length (mm)',side=2,las=0,line=1.4,cex=0.75)
    mtext('Month',side=1,line=1,cex=0.75)
    legend(x=as.Date('2014-06-15'),y=28,c('Season','Constant'),
           lty=1,col=c('black','darkgray'),lwd=1,cex=0.75,bty='n')
    arrows(x0=dateplot.f13[which(L.f13S>=40)[1]],y0=22.5,
           x1=dateplot.f13[which(L.f13S>=40)[1]],y1=40,length=0,lty=2)
    arrows(x0=dateplot.f13[which(L.f13C>=40)[1]],y0=22.5,
           x1=dateplot.f13[which(L.f13C>=40)[1]],y1=40,length=0,lty=2,col='darkgray')
    arrows(x0=dateplot.f13[which(L.f13S>=45)[1]],y0=22.5,
           x1=dateplot.f13[which(L.f13S>=45)[1]],y1=45,length=0,lty=3)
    arrows(x0=dateplot.f13[which(L.f13C>=45)[1]],y0=22.5,
           x1=dateplot.f13[which(L.f13C>=45)[1]],y1=45,length=0,lty=3,col='darkgray')

  #Calculate time difference to 90% aymptotic size
  seas.est <- dateplot.f13[which(L.f13S>=45)[1]]
  const.est <- dateplot.f13[which(L.f13C>=45)[1]]
  (diff.est.d <- as.numeric(seas.est - const.est))
  (diff.est.w <- diff.est.d/7)
  (diff.est.m <- diff.est.d/30)
  #Percent difference:
  (daysto45.s <- as.numeric(dateplot.f13[which(L.f13S>=45)[1]] - dateplot.f13[1]))
  (daysto45.c <- as.numeric(dateplot.f13[which(L.f13C>=45)[1]] - dateplot.f13[1]))
  (daysto45.s - daysto45.c)/daysto45.s*100 
  
  #Calculate time difference to sexuality maturity (40 mm)
  seas.est40 <- dateplot.f13[which(L.f13S>=40)[1]]
  const.est40 <- dateplot.f13[which(L.f13C>=40)[1]]
  (diff.est40.d <- as.numeric(const.est40 - seas.est40))
  (diff.est40.w <- diff.est40.d/7)
  (diff.est40.m <- diff.est40.d/30)
  #Percent difference:
  (daysto40.s <- as.numeric(dateplot.f13[which(L.f13S>=40)[1]] - dateplot.f13[1]))
  (daysto40.c <- as.numeric(dateplot.f13[which(L.f13C>=40)[1]] - dateplot.f13[1]))
  (daysto40.c - daysto40.s)/daysto40.s*100
  
