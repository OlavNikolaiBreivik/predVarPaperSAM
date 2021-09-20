library(stockassessment)

load("code/data/NEAhaddock/model.RData")
dat = fit$data
dat$logobs[dat$aux[,2]==1] = dat$logobs[dat$aux[,2]==1]-3 #make simulation experiment more robust when increasing catch unit
conf = loadConf(dat = dat,"code/conf/NEAhaddock/confHaddock.cfg")
par<-defpar(dat,conf)
fit2 <- sam.fit(dat,conf,par)

conf = fit2$conf
conf$predVarObsLink = matrix(-1, 5, 11)
par<-defpar(dat,conf)
fit <- sam.fit(dat,conf,par)

conf3 = loadConf(dat,file = "code/conf/NEAhaddock/modelFree.cfg")
par3<-defpar(dat,conf3)
fit3 <- sam.fit(dat,conf3,par3)

set.seed(123)
nsim =200
simdataAll <- simulate(fit2, nsim = nsim, full.data = TRUE)
if(sum(is.na(fit2$data$logobs))>0){
  simdata = lapply(simdataAll,function(f){
    f$logobs[which(is.na(fit2$data$logobs))] = NA
    f
  })
}else{
  simdata = simdataAll
}


library(parallel)
ncores = 5
source("code/utils/methods.R")
sim = simStudyPaper(simdata = simdata,fit = fit,nsim = nsim, ncores = ncores, newtonsteps=0)
sim2 = simStudyPaper(simdata = simdata,fit = fit2,nsim = nsim, ncores = ncores, newtonsteps=0)
sim3 = simStudyPaper(simdata = simdata,fit = fit3,nsim = nsim, ncores = ncores, newtonsteps=0)

attr(sim, "fit") <- fit2
attr(sim2, "fit") <- fit2
attr(sim3, "fit") <- fit2


rlFit = as.list(fit2$sdrep,"Est",report = TRUE)
nnn = nsim *length(rlFit$logssb)

meanAbs = rep(NA,nnn)
meanAbs2 = rep(NA,nnn)
meanAbs3 = rep(NA,nnn)

meanAbsF = rep(NA,nnn)
meanAbsF2 = rep(NA,nnn)
meanAbsF3 = rep(NA,nnn)

meanSigma = rep(NA,nnn)
meanSigma2 = rep(NA,nnn)
meanSigma3 = rep(NA,nnn)
meanQL = rep(NA,nnn)
meanQL2 = rep(NA,nnn)
meanQL3 = rep(NA,nnn)
meanQU = rep(NA,nnn)
meanQU2 = rep(NA,nnn)
meanQU3 = rep(NA,nnn)
meanSigmaF = rep(NA,nnn)
meanSigmaF2 = rep(NA,nnn)
meanSigmaF3 = rep(NA,nnn)
meanFQL = rep(NA,nnn)
meanFQL2 = rep(NA,nnn)
meanFQL3 = rep(NA,nnn)
meanFQU = rep(NA,nnn)
meanFQU2 = rep(NA,nnn)
meanFQU3 = rep(NA,nnn)

for(i in 1:nsim){
  c1 = (i-1)*length(rlFit$logssb)+1
  c2 = i*length(rlFit$logssb)

  rl = as.list(sim[[i]]$sdrep,"Est",report = TRUE)
  rlSd = as.list(sim[[i]]$sdrep,"Std",report = TRUE)
  meanAbs[c1:c2] = abs(rlFit$logssb-rl$logssb)
  meanAbsF[c1:c2] = abs(rlFit$logfbar-rl$logfbar)
  meanSigma[c1:c2] = rlSd$logssb
  meanSigmaF[c1:c2] = (rlSd$logfbar)
  meanQL[c1:c2] = (rl$logssb - 1.96*rlSd$logssb) < rlFit$logssb
  meanQU[c1:c2] = (rl$logssb + 1.96*rlSd$logssb) > rlFit$logssb
  meanFQL[c1:c2] = (rl$logfbar - 1.96*rlSd$logfbar) < rlFit$logfbar
  meanFQU[c1:c2] = (rl$logfbar + 1.96*rlSd$logfbar) > rlFit$logfbar


  rl = as.list(sim2[[i]]$sdrep,"Est",report = TRUE)
  rlSd = as.list(sim2[[i]]$sdrep,"Std",report = TRUE)
  meanAbs2[c1:c2] = abs(rlFit$logssb-rl$logssb)
  meanAbsF2[c1:c2] = abs(rlFit$logfbar-rl$logfbar)
  meanSigma2[c1:c2] = rlSd$logssb
  meanSigmaF2[c1:c2] = rlSd$logfbar
  meanQL2[c1:c2] = (rl$logssb - 1.96*rlSd$logssb) < rlFit$logssb
  meanQU2[c1:c2] = (rl$logssb + 1.96*rlSd$logssb) > rlFit$logssb
  meanFQL2[c1:c2] = (rl$logfbar - 1.96*rlSd$logfbar) < rlFit$logfbar
  meanFQU2[c1:c2] = (rl$logfbar + 1.96*rlSd$logfbar) > rlFit$logfbar


  rl = as.list(sim3[[i]]$sdrep,"Est",report = TRUE)
  rlSd = as.list(sim3[[i]]$sdrep,"Std",report = TRUE)
  meanAbs3[c1:c2] = abs(rlFit$logssb-rl$logssb)
  meanAbsF3[c1:c2] = abs(rlFit$logfbar-rl$logfbar)
  meanSigma3[c1:c2] = rlSd$logssb
  meanSigmaF3[c1:c2] = rlSd$logfbar
  meanQL3[c1:c2] = (rl$logssb - 1.96*rlSd$logssb) < rlFit$logssb
  meanQU3[c1:c2] = (rl$logssb + 1.96*rlSd$logssb) > rlFit$logssb
  meanFQL3[c1:c2] = (rl$logfbar - 1.96*rlSd$logfbar) < rlFit$logfbar
  meanFQU3[c1:c2] = (rl$logfbar + 1.96*rlSd$logfbar) > rlFit$logfbar
}

keep = rep(TRUE,nnn)
keep[is.na(meanQL) | is.na(meanQL2) | is.na(meanQL3) ] = FALSE #Not reported confidence intervals
mean(keep) #Proportion with uncertainty intervals


mean(meanAbs[keep])
mean(meanAbs2[keep])
mean(meanAbs3[keep])

mean(meanAbsF[keep])
mean(meanAbsF2[keep])
mean(meanAbsF3[keep])

mean(meanSigma[keep])
mean(meanSigma2[keep])
mean(meanSigma3[keep])
mean(meanSigmaF[keep])
mean(meanSigmaF2[keep])
mean(meanSigmaF3[keep])

mean(meanQL[keep] &meanQU[keep])
mean(meanQL2[keep] &meanQU2[keep])
mean(meanQL3[keep] &meanQU3[keep])

mean(meanFQL[keep] &meanFQU[keep])
mean(meanFQL2[keep] &meanFQU2[keep])
mean(meanFQL3[keep] &meanFQU3[keep])


source("code/utils/plot.R")

par(oma=c(2,3,1,1))
ssbplot(fit2,main ="Haddock SSB",cex.main = 2.5,cex.lab = 1.5,cex.axis=1.3,ylabInternal = "",ci = TRUE)
for(i in 1:nsim){
  ssbplot(sim2[[i]],ci = FALSE, add = TRUE, col = 'green',lwd = 0.5)
}
ssbplot(fit2, ci = FALSE,add = TRUE)
mtext(text="SSB",cex=1.5,side=2,line=0,outer=TRUE)

#Boxplot
trueBeta = exp(fit2$pl$predVarObs)
trueAlpha = exp(fit2$pl$logSdLogObs)
estBeta = matrix(NA,nsim,length(exp(fit2$pl$predVarObs)))
estAlpha = matrix(NA,nsim,length(exp(fit2$pl$predVarObs)))
for(i in 1:nsim){
  estBeta[i,] = exp(sim2[[i]]$pl$predVarObs)
  estAlpha[i,] = exp(sim2[[i]]$pl$logSdLogObs)
}
boxplot(estBeta+1,xlab = "Beta-parameter",ylab = "Estimate",
        main = "Simulation experiment, haddock, beta",cex.main = 2,
        cex.lab = 1.5,
        xaxt = 'n',cex.axis = 1.4)
names = c(expression(beta^("c")),expression(beta[4]^("c")),expression(beta[5-13]^("c")),
          expression(beta^("s1")),expression(beta^("s2")),expression(beta^("s3")),expression(beta^("s4")))
axis(side = 1, at=1:length(names), labels=names,cex.axis = 1.4)
points(trueBeta+1,col = 'red',cex = 2,pty = 3)

boxplot(estAlpha,xlab = "Alpha-parameter",ylab = "Estimate",
        main = "Simulation experiment, haddock, alpha",ylim = c(0,10),cex.main = 2,
        cex.lab = 1.5,
        xaxt = 'n',cex.axis = 1.4)
names = c(expression(alpha^("c")),expression(alpha[4]^("c")),expression(alpha[5-13]^("c")),
          expression(alpha^("s1")),expression(alpha^("s2")),expression(alpha^("s3")),expression(alpha^("s4")))
axis(side = 1, at=1:length(names), labels=names,cex.axis = 1.4)
points(trueAlpha,col = 'red',cex = 2,pty = 3)



