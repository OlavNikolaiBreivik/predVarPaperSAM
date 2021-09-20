library(stockassessment)

load("code/data/NEAcod/model.RData")
dat = fit$data
conf = loadConf(dat = dat,"code/conf/NEAcod/confCod.cfg")
par<-defpar(dat,conf)
fit = sam.fit(dat,conf,par)

conf<-loadConf(dat,"code/conf/NEAcod/modelLink.cfg", patch=TRUE)
par<-defpar(dat,conf)
fit2<-sam.fit(dat,conf,par)

conf<-loadConf(dat,"code/conf/NEAcod/modelFree.cfg", patch=TRUE)
par<-defpar(dat,conf)
fit3<-sam.fit(dat,conf,par)

#aic
AIC(fit,fit2,fit3)

#Parameter estimates
partable(fit2)

#jit
dat = fit$data
dat$logobs = dat$logobs-3 #more robust when using larger catch unit
conf<-loadConf(dat,"code/conf/NEAcod/modelLink.cfg", patch=TRUE)
par<-defpar(dat,conf)
fit2jit <- sam.fit(dat,conf,par)
set.seed(123)
jj = jit(fit2jit,nojit = 50)


#OSA
res = residuals(fit)
res2 = residuals(fit2)
res3 = residuals(fit3)


#Illustrat that stuctures in residuals can be removed
conf<-loadConf(dat,"code/conf/NEAcod/modelLink.cfg", patch=TRUE)
conf$corFlag = 2
par<-defpar(dat,conf)
fit2AR1<-sam.fit(dat,conf,par)
res2AR1 = residuals(fit2AR1)


#See extra script for simstudy

#Plots in paper
source("code/utils/plot.R") #Minor change in plotting functionality

par(oma=c(2,2,1,1))
ssbplot(fit2,main ="Cod SSB",cex.main = 2.5,cex.lab = 1.5,cex.axis=1.3,ylabInternal = "",lwd = 7)
ssb1 = exp(as.list(fit$sdrep,"Est",report = TRUE)$logssb)
lines(rev(fit$data$years), rev(ssb1), lwd=7, lty=2, col="red")
ssb3 = exp(as.list(fit3$sdrep,"Est",report = TRUE)$logssb)
lines(rev(fit$data$years), rev(ssb3), lwd=7, lty=3, col="blue")
mtext(text="SSB",cex=1.5,side=2,line=0.2,outer=TRUE)


par(oma=c(1,1.5,1,1))
fbarplot(fit2,partial = FALSE,main ="Cod fishing mortality",cex.main = 2.5,cex.lab = 1.3,cex.axis=1.5,ylabInternal = "",lwd = 7)
fbar1 = exp(as.list(fit$sdrep,"Est",report = TRUE)$logfbar)
fbar1 = fbar1[-length(fbar1)]
years = fit$data$years
years = years[-length(years)]
lines(rev(years), rev(fbar1), lwd=7, lty=2, col="red")
fbar3 = exp(as.list(fit3$sdrep,"Est",report = TRUE)$logfbar)
fbar3 = fbar3[-length(fbar3)]
lines(rev(years), rev(fbar3), lwd=7, lty=3, col="blue")
fbarlab=substitute(bar(F)[X-Y],list(X=5,Y=10))
mtext(text=fbarlab,cex=1.5,side=2,line=-1,outer=TRUE)


par(oma=c(1,3,1,1))
catchplot(fit,main = "Aggregated catch cod model 1",cex.main = 2,cex.lab = 1.5,cex.axis=1.5,ylabInternal = "",ylim = c(0,1500000))
mtext(text="Catch (tons)",cex=1.5,side=2,line=1.3,outer=TRUE)

par(oma=c(1,3,1,1))
catchplot(fit2,main = "Aggregated catch cod model 2",cex.main = 2,cex.lab = 1.5,cex.axis=1.5,ylabInternal = "",ylim = c(0,1500000))
mtext(text="Catch (tons)",cex=1.5,side=2,line=1.3,outer=TRUE)

par(oma=c(1,3,1,1))
catchplot(fit3,main = "Aggregated catch cod model 3",cex.main = 2,cex.lab = 1.5,cex.axis=1.5,ylabInternal = "",ylim = c(0,1500000))
mtext(text="Catch (tons)",cex=1.5,side=2,line=1.3,outer=TRUE)




plotRes(x = res,zmax = 2,axes = FALSE,cex.axis = 1.5)
title("OSA residuals cod model 1",line = 1,cex.main=2.2)
mtext(text="Age",cex=2,side=2,line=2.2,outer=FALSE)
mtext(text="Year",cex=2,side=1,line=2.2,outer=FALSE)
plotRes(x = res2,zmax = 2,axes = FALSE,cex.axis = 1.5)
title("OSA residuals cod model 2",line = 1,cex.main=2.2)
mtext(text="Age",cex=2,side=2,line=2,outer=FALSE)
mtext(text="Year",cex=2,side=1,line=2,outer=FALSE)
plotRes(x = res3,zmax = 2,axes = FALSE,cex.axis = 1.5)
title("OSA residuals cod model 3",line = 1,cex.main=2.2)
mtext(text="Age",cex=2,side=2,line=2,outer=FALSE)
mtext(text="Year",cex=2,side=1,line=2,outer=FALSE)
plotRes(x = res2AR1,zmax = 2,axes = FALSE,cex.axis = 1.5)
title("OSA residuals cod adjusted model 2",line = 1,cex.main=2.2)
mtext(text="Age",cex=2,side=2,line=2,outer=FALSE)
mtext(text="Year",cex=2,side=1,line=2,outer=FALSE)


#Trompetplot
hh = which(res$fleet ==1 & res$age>4)
plot(res$mean[hh] ,res$residual[hh],ylim = c(-4,4.5),xlim = c(-2,13),
     ylab = "OSA residual", xlab = "Predicted log-observation",
     main = "Catch of age 5+ cod, model 1",
     cex.main = 1.8,cex.lab = 1.5,cex = 1.4,cex.axis = 1.4)
abline(h = 0)
hh = which(res2$fleet ==1 & res2$age>4)
plot(res2$mean[hh] ,res2$residual[hh],ylim = c(-4,4.5),xlim = c(-2,13),
     ylab = "OSA residual", xlab = "Predicted log-observation",
     main = "Catch of age 5+ cod, model 2",
     cex.main = 1.8,cex.lab = 1.5,cex = 1.4,cex.axis = 1.4)
abline(h = 0)




plot(res,zmax = 3)
mtext("OSA residuals model 1 cod", line = 2,at = 0, cex = 2.5)
plot(res2,zmax = 3)
mtext("OSA residuals model 2 cod", line = 2,at = 0, cex = 2.5)
plot(res3,zmax = 3)
mtext("OSA residuals model 3 cod", line = 2,at = 0, cex = 2.5)


#Sd vs pred, currently not in article.
N = 20000
ss = rmvnorm(N,mu = fit2$opt$par,Sigma = fit2$sdrep$cov.fixed)
alpha = exp(ss[,which(names(fit2$opt$par)=="logSdLogObs")])
beta = exp(ss[,which(names(fit2$opt$par)=="predVarObs")])

plotLink = function(fitSandard,fitLink,f,a,col = 1,main = "",ylim = c(0,3)){

  aaConstant = fitSandard$conf$keyVarObs[f,a]+1

  aaa = fitLink$conf$keyVarObs[f,a]+1
  aab = fitLink$conf$predVarObs[f,a]+1

  age = which(fitSandard$conf$keyVarObs[f,]==aaConstant-1) + fitSandard$conf$minAge -1
  obsIndex = which(fitLink$data$aux[,2]==f & fitLink$data$aux[,3] %in% age)

  pred = fitLink$rep$predObs
  pred = pred[obsIndex]
  obs = fitLink$data$logobs[obsIndex]
  obs = obs[!is.na(obs)]

  l = seq(min(pred),max(pred), by = 0.01)

  sd = matrix(0,N,length(l))
  for(i in 1:N){
    if(aab!=0){
      sd[i,] = sqrt(log(alpha[i,aaa]* exp(l*(beta[i,aab]-1)) +1))
    }else{
      sd[i,] = alpha[i,aaa]
    }
  }

  sdM = rep(0,length(l))
  sdL = rep(0,length(l))
  sdU = rep(0,length(l))
  for(i in 1:length(l)){
    sdM[i] = mean(sd[,i])
    sdL[i] = quantile(sd[,i], 0.025)
    sdU[i] = quantile(sd[,i], 0.975)
  }

  if(length(age)==1){
    plot(l,sdM,type = 'l', ylim = ylim,lwd = 2,col = col,main = paste0(main, " age ",age),
         ylab = "", xlab = "",cex.main = 2.5,
         cex.lab = 1.5,
         cex.axis = 1.5)
  }else{
    plot(l,sdM,type = 'l', ylim = ylim,lwd = 2,col = col,main = paste0(main, " age ",min(age),"-",max(age)),
         ylab = "", xlab = "",cex.main = 2.5,
         cex.lab = 2,
         cex.axis = 1.5)
  }
  mtext("Standard deviation", cex = 2,side = 2, line=2.7,)
  mtext("Log-prediction", cex = 2,side = 1, line=2.7,)
  lines(l,sdL,col = col)
  lines(l,sdU,col = col)
  abline(h = exp(fitSandard$pl$logSdLogObs[aaConstant]), col = col+1)
  abline(h = exp(fitSandard$pl$logSdLogObs[aaConstant] + 1.96*fitSandard$plsd$logSdLogObs[aaConstant]),col = col+1,lt = 2)
  abline(h = exp(fitSandard$pl$logSdLogObs[aaConstant] - 1.96*fitSandard$plsd$logSdLogObs[aaConstant]),col = col+1,lt = 2)

  points(pred,rep(exp(fitSandard$pl$logSdLogObs[aaConstant]),length(pred)))
}
ylim = c(0,1.3)
plotLink(fitSandard = fit,fitLink = fit2,f = 1,a = 1,main = "Catch sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 1,a = 2,main = "Catch sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 1,a = 3,main = "Catch sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 1,a = 9,main = "Catch sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 1,a = 11,main = "Catch sd",ylim = ylim)

plotLink(fitSandard = fit,fitLink = fit2,f = 2,a = 1,main = "Survey 1 sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 2,a = 2,main = "Survey 1 sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 2,a = 6,main = "Survey 1 sd",ylim = ylim)

plotLink(fitSandard = fit,fitLink = fit2,f = 3,a = 1,main = "Survey 2 sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 3,a = 2,main = "Survey 2 sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 3,a = 6,main = "Survey 2 sd",ylim = ylim)

plotLink(fitSandard = fit,fitLink = fit2,f = 4,a = 1,main = "Survey 3 sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 4,a = 7,main = "Survey 3 sd",ylim = ylim)

plotLink(fitSandard = fit,fitLink = fit2,f = 6,a = 1,main = "Survey 4 sd",ylim = ylim)

