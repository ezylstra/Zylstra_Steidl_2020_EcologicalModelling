##############################################################################################
#Zylstra, ER, and RJ Steidl. 2020. A Bayesian state-space model for seasonal growth of 
#terrestrial vertebrates. Ecological Modelling 420:108975

#Estimating growth of canyon treefrogs using capture-recapture data with a 
#seasonal version of a von Bertalanffy, length-at-first-capture (LFC) growth model.
##############################################################################################

#------------------------------------------------------------------------------------------# 
# Note that notation below is different than that used in the manuscript.
# Here, asymptotic length is represented by Linf (in the manuscript, we used L)
# and length at first capture is represented by L1 (in the manuscript, we use lambda)
#------------------------------------------------------------------------------------------#


#Set working directory...

# install.packages(c('plyr','jagsUI'))
library(plyr)
library(jagsUI)

#------------------------------------------------------------#
# Read in data and format date
#------------------------------------------------------------#
  ctf <- read.table("AllFrogs.csv",sep=",",header=TRUE,na.strings=c("NA",""),strip.white=T,stringsAsFactors=F)
  
  ctf$CapDate <- as.Date(ctf$CapDate,format="%Y-%m-%d")
  ctf$Year <- as.numeric(format(ctf$CapDate,'%Y'))

#------------------------------------------------------------#
# Prepping data for growth analysis
#------------------------------------------------------------#
#Remove data from individuals that were not marked (in 2016)
  ctfsub <- ctf[!is.na(ctf$Tag),]

#Remove individuals were were never recaptured
  ncaps <- ddply(ctfsub,.(Tag),summarize,ncap=length(CapDate))
  recaps <- ncaps$Tag[ncaps$ncap>1]
  ctfsub <- ctfsub[ctfsub$Tag %in% recaps,]

#Create a dataframe, where each row represents one frog
  df.frog <- ddply(ctfsub,.(Tag),summarize,ncap=length(CapDate),Sex=Sex[1],CapDate1=min(CapDate),Year1=min(Year))
  df.frog$Fem <- ifelse(df.frog$Sex=='f',1,ifelse(df.frog$Sex=='m',0,NA))
  #Calculate the time of first capture (in fractions of a year)
  df.frog$t <- as.numeric(df.frog$CapDate1 - as.Date(paste(df.frog$Year1,'01','01',sep='-'),format='%Y-%m-%d'))/365
  #Create an index for individuals
  df.frog <- df.frog[with(df.frog,order(Fem,na.last=F,CapDate1)),]
  df.frog$ID <- 1:nrow(df.frog)
  
#Create a dataframe, where each row represents a capture
  df.cap <- ctfsub[,c('Tag','CapDate','SVL','Year')]
  df.cap <- join(df.cap,df.frog[,c('Tag','Fem','CapDate1','Year1','t','ID')],by='Tag',type='left')
  df.cap <- df.cap[with(df.cap,order(ID,CapDate)),]
  #Calculate the time elapsed since first capture (in fractions of a year)
  df.cap$d <- as.numeric(df.cap$CapDate - df.cap$CapDate1)/365

  nind <- nrow(df.frog)
  nunk <- sum(is.na(df.frog$Fem))
  nmeas <- nrow(df.cap)
  
#-------------------------------------------------------------------------------------#
# JAGS model
# Seasonal LFC model: Estimate C, no individual heterogeneity in k, sex-specific SD for Linf, posterior predictive checks
#-------------------------------------------------------------------------------------#
  sink('CTFgrowth_LFC_Season_EstC.txt')
  cat("
    model {

    for(m in 1:nmeas){
      y[m] ~ dnorm(mu[m],tau.meas)
      mu[m] <- L1[ind[m]] + (Linf[ind[m]]-L1[ind[m]])*(1-exp(-k[ind[m]]*d[m] -
               C*k[ind[m]]/(2*pi)*sin(2*pi*(t[m]+d[m]-ts)) +
               C*k[ind[m]]/(2*pi)*sin(2*pi*(t[m]-ts))))
      y.pred[m] ~ dnorm(mu[m],tau.meas)
    }

    for(i in 1:nind){
      L1[i] ~ dnorm(alpha[1] + alpha[2]*fem[i],tau.alpha)T(0,Linf[i])
      Linf[i] ~ dnorm(beta[1] + beta[2]*fem[i],tau.beta[fem[i]+1])T(0,)
      log(k[i]) <- gamma[1] + gamma[2]*fem[i]
      fem[i] ~ dbern(psi)
    }

    alpha[1] ~ dnorm(0,0.001)
    alpha[2] ~ dnorm(0,0.001)
    beta[1] ~ dnorm(0,0.001)
    beta[2] ~ dnorm(0,0.001)
    gamma[1] ~ dnorm(0,0.001)
    gamma[2] ~ dnorm(0,0.001)

    tau.meas <- 1/(sd.meas*sd.meas)
    sd.meas ~ dt(0,pow(5,-2),1)T(0,)

    tau.alpha <- 1/(sd.alpha*sd.alpha)
    sd.alpha ~ dt(0,pow(5,-2),1)T(0,)

    for(g in 1:2){
      tau.beta[g] <- 1/(sd.beta[g]*sd.beta[g])
      sd.beta[g] ~ dt(0,pow(5,-2),1)T(0,)
    }

    psi ~ dunif(0,1)

    pi <- 3.14159265359
    C ~ dunif(0,1)       #Setting C = 1 assumes seasonal cessation of growth
    ts ~ dunif(0,1)

    Linf.fem <- beta[1] + beta[2]
    L1.fem <- alpha[1] + alpha[2]
    k.male <- exp(gamma[1])
    k.fem <- exp(gamma[1] + gamma[2])
    ts.julian <- ts*365

  } #model
  ", fill=TRUE)
  sink()

#------------------------------------------------------------#
# Bundle items for JAGS
#------------------------------------------------------------#
  jdata <- list(y=df.cap$SVL,
                d=df.cap$d,
                t=df.cap$t,
                ind=df.cap$ID,
                fem=df.frog$Fem,
                nind=nind,
                nmeas=nmeas)

  inits <- function(){list(alpha=c(runif(1,30,60),runif(1,-2,2)),
                           beta=c(runif(1,40,60),runif(1,-2,2)),
                           gamma=c(runif(1,-2,2),runif(1,-2,2)),
                           sd.meas=runif(1,0,3),
                           sd.alpha=runif(1,0,5),
                           sd.beta=runif(2,0,5),
                           C=runif(1,0,1),
                           ts=runif(1,0,1),
                           psi=dunif(1,0,1),
                           fem=c(rep(1,nunk),rep(NA,nind-nunk)))}

  params <- c('alpha','beta','C','ts','ts.julian','L1.fem','Linf.fem','k.male','k.fem',
              'sd.meas','sd.alpha','sd.beta','psi','y.pred')

  nc <- 3
  na <- 20000  
  nb <- 80000  
  ni <- 40000  
  nt <- 20     
  ni.tot <- nb + ni
  
  jags.model <- 'CTFgrowth_LFC_Season_EstC.txt'
  
#------------------------------------------------------------#
# Run model and look at results, assess fit
#------------------------------------------------------------#
  fit <- jags(data=jdata, inits=inits, parameters.to.save=params, model.file=jags.model,
              n.chains=nc, n.adapt=na, n.iter=ni.tot, n.burnin=nb, n.thin=nt,
              parallel=T, store.data=T)
  
  print(fit)
  
  out <- fit$samples
  plot(out,ask=T)
  
#Percent of samples from posterior distribution of C that are near 1
  C <- unlist(out[,'C'])
  sum(C>=0.90)/length(C)
  sum(C>=0.95)/length(C)
  
#Evaluation of fit 
  outm <- do.call(rbind,out)
  y.ind <- grep('y.pred',colnames(out[[1]]))
  outm.ppc <- outm[,y.ind]

  y.mn <- apply(outm.ppc,1,mean)
    p.mn <- sum(y.mn>=mean(df.cap$svl))/length(y.mn)
  y.mx <- apply(outm.ppc,1,max)
    p.mx <- sum(y.mx>=max(df.cap$svl))/length(y.mx)
  y.sd <- apply(outm.ppc,1,sd)
    p.sd <- sum(y.sd>=sd(df.cap$svl))/length(y.sd)
    
#Three-panel figure with posterior predictive distributions:
  par(mfrow=c(1,3),mar=c(3,0.5,0.5,0),oma=c(0,0,0,0.5),mgp=c(3,0.5,0),tcl=-0.25)  
    hist(y.mn,breaks=16,col='gray',main='',axes=F,xlim=c(43.25,43.75))  
      axis(side=1)
      mtext('Mean',side=1,line=2)
      abline(v=mean(df.cap$svl),col='blue',lwd=lwidth) 
    hist(y.mx,breaks=15,col='gray',main='',axes=F,xlim=c(51,60))  
      axis(side=1)
      mtext('Maximum',side=1,line=2)
      abline(v=max(df.cap$svl),col='blue',lwd=lwidth) 
    hist(y.sd,breaks=14,col='gray',main='',axes=F,xlim=c(3.4,3.9))  
      axis(side=1)
      mtext('SD',side=1,line=2)
      abline(v=sd(df.cap$svl),col='blue',lwd=lwidth) 

#------------------------------------------------------------------------------#
# Create figure with predicted growth curves (and histograms with capture data)
#------------------------------------------------------------------------------#
#First, prep the length-frequncy data from ALL captures 
#(ctf is the original dataframe that includes individuals that were not marked or recaptured)
  
  #Survey dates (approximate midpoints)
  sdates <- as.Date(c('2014-04-15','2014-05-15','2014-06-15','2014-09-30','2014-10-30',
                      '2015-04-15','2015-05-15','2015-06-15','2015-09-30','2015-10-30','2016-06-05')) 
  #Getting data from histograms of SVL by occasion for each sex (histogram objects named f1, f2, m1, m2, etc.)
  #Removing few observations with SVL>54 for figure
  #Using Sex40 so that all individuals <40 mm are listed as unknown, even if sex was determined with subsequent capture when SVL >= 40
  for(i in 1:11){
    assign(paste('f',i,sep=''),hist(ctf$SVL[which(ctf$Occasion==i & ctf$SVL>=40 & ctf$SVL<55 & ctf$Sex40=='f')],breaks=seq(34.5,54.5,by=1)))
    assign(paste('m',i,sep=''),hist(ctf$SVL[which(ctf$Occasion==i & ctf$SVL>=40 & ctf$SVL<55 & ctf$Sex40=='m')],breaks=seq(34.5,54.5,by=1)))
    assign(paste('u',i,sep=''),hist(ctf$SVL[which(ctf$Occasion==i & ctf$SVL>34 & ctf$SVL<40 & ctf$Sex40=='u')],breaks=seq(34.5,54.5,by=1)))
  }
  
  breaks <- f1$breaks
  fcounts <- matrix(c(f1$counts,f2$counts,f3$counts,f4$counts,f5$counts,
                      f6$counts,f7$counts,f8$counts,f9$counts,f10$counts,f11$counts),
                    byrow=T,nrow=11)
  mcounts <- matrix(c(m1$counts,m2$counts,m3$counts,m4$counts,m5$counts,
                      m6$counts,m7$counts,m8$counts,m9$counts,m10$counts,m11$counts),
                    byrow=T,nrow=11)
  ucounts <- matrix(c(u1$counts,u2$counts,u3$counts,u4$counts,u5$counts,
                      u6$counts,u7$counts,u8$counts,u9$counts,u10$counts,u11$counts),
                    byrow=T,nrow=11)
  ucounts <- round(ucounts/2)

#Extract mean of posterior distributions for model parameters and assign other scalar values
  Linf.mS <- mean(outm[,'beta[1]'])
  Linf.fS <- mean(outm[,'Linf.fem'])
  k.mS <- mean(outm[,'k.male'])
  k.fS <- mean(outm[,'k.fem'])
  ts <- mean(outm[,'ts'])
  C <- 1
  Lt1 <- 22.5  #size at metamorphosis
  int <- 1/365

#Create vectors to fill with predicted lengths
  yrdays <- c(rep(1:365,3),1:183)  #2013 through 1 July 2016 
  #Setting metamorphosis dates:
  stdate.f13 <- as.Date('2013-09-01'); stdate.f14 <- as.Date('2014-09-01'); stdate.f15 <- as.Date('2015-09-01')
  stdate.s13 <- as.Date('2013-07-01'); stdate.s14 <- as.Date('2014-07-01'); stdate.s15 <- as.Date('2015-07-01')
  enddate <- as.Date('2016-07-01')
  #converting metamorphosis dates to fraction of a year
  day1.f13 <- as.numeric(strftime(stdate.f13,format='%j')); t1.f13 <- day1.f13/365
  day1.f14 <- as.numeric(strftime(stdate.f14,format='%j')); t1.f14 <- day1.f14/365
  day1.f15 <- as.numeric(strftime(stdate.f15,format='%j')); t1.f15 <- day1.f15/365
  day1.s13 <- as.numeric(strftime(stdate.s13,format='%j')); t1.s13 <- day1.s13/365
  day1.s14 <- as.numeric(strftime(stdate.s14,format='%j')); t1.s14 <- day1.s14/365
  day1.s15 <- as.numeric(strftime(stdate.s15,format='%j')); t1.s15 <- day1.s15/365
  #create vectors for length of individuals in each cohort over time
  l13 <- length(yrdays); l14 <- l13-365; l15 <- l13-730
  Lf.f13S <- c(Lt1,rep(NA,l13-day1.f13))
  Lf.f14S <- c(Lt1,rep(NA,l14-day1.f14))
  Lf.f15S <- c(Lt1,rep(NA,l15-day1.f15))
  Lf.s13S <- c(Lt1,rep(NA,l13-day1.s13))
  Lf.s14S <- c(Lt1,rep(NA,l14-day1.s14))
  Lf.s15S <- c(Lt1,rep(NA,l15-day1.s15))
  Lm.f13S <- c(Lt1,rep(NA,l13-day1.f13))
  Lm.f14S <- c(Lt1,rep(NA,l14-day1.f14))
  Lm.f15S <- c(Lt1,rep(NA,l15-day1.f15))
  Lm.s13S <- c(Lt1,rep(NA,l13-day1.s13))
  Lm.s14S <- c(Lt1,rep(NA,l14-day1.s14))
  Lm.s15S <- c(Lt1,rep(NA,l15-day1.s15))

#Predict length for each day after metamorphosis based on model estimates:
  for(i in 2:(l13-day1.f13)){
    Lf.f13S[i] <- Lt1 + (Linf.fS-Lt1)*(1-exp(-k.fS*int*(i-1) - C*k.fS/(2*pi)*sin(2*pi*(t1.f13+int*(i-1)-ts)) +
                  C*k.fS/(2*pi)*sin(2*pi*(t1.f13-ts)))) 
    Lm.f13S[i] <- Lt1 + (Linf.mS-Lt1)*(1-exp(-k.mS*int*(i-1) - C*k.mS/(2*pi)*sin(2*pi*(t1.f13+int*(i-1)-ts)) +
                  C*k.mS/(2*pi)*sin(2*pi*(t1.f13-ts))))
  }
  for(i in 2:(l14-day1.f14)){
    Lf.f14S[i] <- Lt1 + (Linf.fS-Lt1)*(1-exp(-k.fS*int*(i-1) - C*k.fS/(2*pi)*sin(2*pi*(t1.f14+int*(i-1)-ts)) +
                  C*k.fS/(2*pi)*sin(2*pi*(t1.f14-ts)))) 
    Lm.f14S[i] <- Lt1 + (Linf.mS-Lt1)*(1-exp(-k.mS*int*(i-1) - C*k.mS/(2*pi)*sin(2*pi*(t1.f14+int*(i-1)-ts)) +
                  C*k.mS/(2*pi)*sin(2*pi*(t1.f14-ts)))) 
  }
  for(i in 2:(l15-day1.f15)){
    Lf.f15S[i] <- Lt1 + (Linf.fS-Lt1)*(1-exp(-k.fS*int*(i-1) - C*k.fS/(2*pi)*sin(2*pi*(t1.f15+int*(i-1)-ts)) +
                  C*k.fS/(2*pi)*sin(2*pi*(t1.f15-ts)))) 
    Lm.f15S[i] <- Lt1 + (Linf.mS-Lt1)*(1-exp(-k.mS*int*(i-1) - C*k.mS/(2*pi)*sin(2*pi*(t1.f15+int*(i-1)-ts)) +
                  C*k.mS/(2*pi)*sin(2*pi*(t1.f15-ts)))) 
  }
  for(i in 2:(l13-day1.s13)){  
    Lf.s13S[i] <- Lt1 + (Linf.fS-Lt1)*(1-exp(-k.fS*int*(i-1) - C*k.fS/(2*pi)*sin(2*pi*(t1.s13+int*(i-1)-ts)) +
                  C*k.fS/(2*pi)*sin(2*pi*(t1.s13-ts)))) 
    Lm.s13S[i] <- Lt1 + (Linf.mS-Lt1)*(1-exp(-k.mS*int*(i-1) - C*k.fS/(2*pi)*sin(2*pi*(t1.s13+int*(i-1)-ts)) +
                  C*k.mS/(2*pi)*sin(2*pi*(t1.s13-ts)))) 
  }
  for(i in 2:(l14-day1.s14)){  
    Lf.s14S[i] <- Lt1 + (Linf.fS-Lt1)*(1-exp(-k.fS*int*(i-1) - C*k.fS/(2*pi)*sin(2*pi*(t1.s14+int*(i-1)-ts)) +
                  C*k.fS/(2*pi)*sin(2*pi*(t1.s14-ts)))) 
    Lm.s14S[i] <- Lt1 + (Linf.mS-Lt1)*(1-exp(-k.mS*int*(i-1) - C*k.fS/(2*pi)*sin(2*pi*(t1.s14+int*(i-1)-ts)) +
                  C*k.mS/(2*pi)*sin(2*pi*(t1.s14-ts)))) 
  }
  for(i in 2:(l15-day1.s15)){  
    Lf.s15S[i] <- Lt1 + (Linf.fS-Lt1)*(1-exp(-k.fS*int*(i-1) - C*k.fS/(2*pi)*sin(2*pi*(t1.s15+int*(i-1)-ts)) +
                  C*k.fS/(2*pi)*sin(2*pi*(t1.s15-ts)))) 
    Lm.s15S[i] <- Lt1 + (Linf.mS-Lt1)*(1-exp(-k.mS*int*(i-1) - C*k.fS/(2*pi)*sin(2*pi*(t1.s15+int*(i-1)-ts)) +
                  C*k.mS/(2*pi)*sin(2*pi*(t1.s15-ts)))) 
  }
  
#Dates for plotting
  dateplot.f13 <- seq(from=stdate.f13,to=enddate,by='days')
  dateplot.f14 <- seq(from=stdate.f14,to=enddate,by='days')
  dateplot.f15 <- seq(from=stdate.f15,to=enddate,by='days')
  dateplot.s13 <- seq(from=stdate.s13,to=enddate,by='days')
  dateplot.s14 <- seq(from=stdate.s14,to=enddate,by='days')
  dateplot.s15 <- seq(from=stdate.s15,to=enddate,by='days')
 
#Stacked plot (Females top, Males bottom) ranging from 1 Apr 2014 - 1 July 2016
  par(mfrow=c(2,1),mar=c(0.5,2.6,0.5,0.6)+0.1,oma=c(1.6,0,0,0))
  plot(1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',
       xlim=c(dateplot.s13[275]-10,enddate),ylim=c(21,54),bty='n',xaxs='i',yaxs='i')  #185 is 1 Jan 2014
    usr <- par('usr')  #these are plotting limits (incl extra bit)
    axis(1,at=c(usr[1],usr[2]),tck=F,labels=F)
    axis.Date(1,at=seq(as.Date('2014-04-01'),enddate,by='quarter'),tcl=-0.2,labels=F)
    axis(2,at=c(usr[3],usr[4]),tck=F,labels=F)
    axis(2,at=seq(22,54,by=4),labels=seq(22,54,by=4),las=1,
         tcl=-0.2,cex.axis=0.9,mgp=c(1.5,0.4,0))
    #Add in gray lines at survey dates
    for(i in 1:length(sdates)){
      arrows(x0=sdates[i],y0=22.5,x1=sdates[i],y1=54,length=0,col='gray')
    }
    #Add in histograms with SVL by occasion (females and unknowns[half])
    for(i in 1:11){
      rect(sdates[i]-fcounts[i,]*0.66,breaks[-length(breaks)],sdates[i],breaks[-1],
           col=rgb(205/255,91/255,69/255,alpha=0.4))
      rect(sdates[i]-ucounts[i,]*0.66,breaks[-length(breaks)],sdates[i],breaks[-1])
    }
    lines(Lf.s13S[275:length(Lf.s13S)]~dateplot.s13[275:length(dateplot.s13)],lty=2,lwd=1)
    lines(Lf.s14S~dateplot.s14,lty=2,lwd=1)
    lines(Lf.s15S~dateplot.s15,lty=2,lwd=1)
    lines(Lf.f13S[213:length(Lf.f13S)]~dateplot.f13[213:length(dateplot.f13)],lwd=1)
    lines(Lf.f14S~dateplot.f14,lwd=1)
    lines(Lf.f15S~dateplot.f15,lwd=1)
    mtext('SVL (mm)',side=2,las=0,line=1.7,cex=0.9)
    legend(x=as.Date('2016-01-05'),y=27,c('Summer','Fall'),
           lty=c(2,1),lwd=1,cex=0.87,bty='n')
    text(x=as.Date('2016-05-25'),y=52.5,'Females',cex=0.93,adj=c(1,0))
  plot(1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',
       xlim=c(dateplot.s13[275]-10,enddate),ylim=c(21,54),bty='n',xaxs='i',yaxs='i')
    usr <- par('usr')  #these are plotting limits (incl extra bit)
    axis(1,at=c(usr[1],usr[2]),tck=F,labels=F)
    axis.Date(1,at=seq(as.Date('2014-04-01'),enddate,by='quarter'),format='%b',
            tcl=-0.2,cex.axis=0.9,mgp=c(1.5,0.1,0))
    axis(2,at=c(usr[3],usr[4]),tck=F,labels=F)
    axis(2,at=seq(22,54,by=4),labels=seq(22,54,by=4),las=1,
         tcl=-0.2,cex.axis=0.9,mgp=c(1.5,0.4,0))
    #Add in gray lines at survey dates
    for(i in 1:length(sdates)){
      arrows(x0=sdates[i],y0=22.5,x1=sdates[i],y1=54,length=0,col='gray')
    }
    #Add in histograms with SVL by occasion (males and unknowns[all])
    for(i in 1:11){
      rect(sdates[i]-mcounts[i,]*0.66,breaks[-length(breaks)],sdates[i],breaks[-1],
           col=rgb(24/255,116/255,205/255,alpha=0.4))  #could use border=F
      rect(sdates[i]-ucounts[i,]*0.66,breaks[-length(breaks)],sdates[i],breaks[-1])
    }
    lines(Lm.s13S[275:length(Lm.s13S)]~dateplot.s13[275:length(dateplot.s13)],lty=2,lwd=1)
    lines(Lm.s14S~dateplot.s14,lty=2,lwd=1)
    lines(Lm.s15S~dateplot.s15,lty=2,lwd=1)
    lines(Lm.f13S[213:length(Lm.f13S)]~dateplot.f13[213:length(dateplot.f13)],lwd=1)
    lines(Lm.f14S~dateplot.f14,lwd=1)
    lines(Lm.f15S~dateplot.f15,lwd=1)
    mtext('2014',side=1,outer=T,line=0.6,at=0.17,adj=0,cex=0.9)
    mtext('2015',side=1,outer=T,line=0.6,at=0.51,adj=0,cex=0.9)
    mtext('2016',side=1,outer=T,line=0.6,at=0.85,adj=0,cex=0.9)
    mtext('SVL (mm)',side=2,las=0,line=1.7,cex=0.9)
    text(x=as.Date('2016-05-25'),y=52.5,'Males',cex=0.93,adj=c(1,0))

#------------------------------------------------------------------------------------------#
# Calculating time required to reach 40 mm (sexual maturity) or 90/95% of asymptotic length
#------------------------------------------------------------------------------------------#
  Linf <- c(Linf.fS,Linf.mS)
  Linf90 <- Linf*0.90
  Linf95 <- Linf*0.95
  sexmat <- c(40,40)
  Linftime <- data.frame(sex=rep(c('f','m'),2),season=c('s','s','f','f'))
  predictlist <- paste('L',Linftime$sex,'.',Linftime$season,'13S',sep='')
  for(i in c(1,3)){
    Linftime$t90[i] <- which(get(predictlist[i])>Linf90[1])[1]
    Linftime$t90[i+1] <- which(get(predictlist[i+1])>Linf90[2])[1]
    Linftime$t95[i] <- which(get(predictlist[i])>Linf95[1])[1]
    Linftime$t95[i+1] <- which(get(predictlist[i+1])>Linf95[2])[1]
    Linftime$t40[i] <- which(get(predictlist[i])>sexmat[1])[1]
    Linftime$t40[i+1] <- which(get(predictlist[i+1])>sexmat[2])[1]
  }
  Linftime$t90.mon <- Linftime$t90/30
  Linftime$t95.mon <- Linftime$t95/30
  Linftime$t40.mon <- Linftime$t40/30
  Linftime$t90.yr <- Linftime$t90/365
  Linftime$t95.yr <- Linftime$t95/365
  Linftime$t40.yr <- Linftime$t40/365
  Linftime
  
  #How long for males from summer cohort to reach sexual maturity?
  Linftime$t40.mon[2]
  #How long for males from fall cohort to reach sexual maturity?
  Linftime$t40.mon[4]  
  #How long for females from summer cohort to reach sexual maturity?
  Linftime$t40.mon[1]
  #How long for females from fall cohort to reach sexual maturity?
  Linftime$t40.mon[3]
  
  #How much longer for males to reach sexual maturity (assuming both sexes mature at 40)?
  Linftime$t40.mon[c(2,4)]-Linftime$t40.mon[c(1,3)]  

