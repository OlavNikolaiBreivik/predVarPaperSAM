library(stockassessment)
load("code/data/NEAhaddock/model.RData")
dat = fit$data
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

#aic
AIC(fit,fit2,fit3)

#jit
set.seed(123)
jj = jit(fit2,nojit = 50,ncores = 5)

#OSA
res = residuals(fit)
res2 = residuals(fit2)
res3 = residuals(fit3)

#sim, see simulationExperiment.R

#Parameter estimates
partable(fit2)

source("code/utils/plot.R")#minor modifications of plotting functions

#Construct plots for haddock in paper
par(oma=c(2,2,1,1))
ssbplot(fit2,main ="Haddock SSB",cex.main = 2.5,cex.lab = 1.5,cex.axis=1.3,ylabInternal = "",lwd = 7)
ssb1 = exp(as.list(fit$sdrep,"Est",report = TRUE)$logssb)
lines(rev(fit$data$years), rev(ssb1), lwd=7, lty=2, col="red")
ssb3 = exp(as.list(fit3$sdrep,"Est",report = TRUE)$logssb)
lines(rev(fit$data$years), rev(ssb3), lwd=7, lty=3, col="blue")
legend(x=1950,y=650000,legend = c("Model 1", "Model 2", "Model 3")
       ,col=c("red","black","blue"),lty = c(2,1,3),lwd = 6,box.lty =0,cex =1.7)
mtext(text="SSB",cex=1.5,side=2,line=0.2,outer=TRUE)



par(oma=c(1,1.5,1,1))
fbarplot(fit2,partial = FALSE,main ="Haddock fishing mortality",cex.main = 2.5,cex.lab = 1.3,cex.axis=1.5,ylabInternal = "",lwd = 7)
fbar1 = exp(as.list(fit$sdrep,"Est",report = TRUE)$logfbar)
fbar1 = fbar1[-length(fbar1)]
years = fit$data$years
years = years[-length(years)]
lines(rev(years), rev(fbar1), lwd=7, lty=2, col="red")
fbar3 = exp(as.list(fit3$sdrep,"Est",report = TRUE)$logfbar)
fbar3 = fbar3[-length(fbar3)]
lines(rev(years), rev(fbar3), lwd=7, lty=3, col="blue")
legend(x=1993,y=0.9,legend = c("Model 1", "Model 2","Model 3")
       ,col=c("red","black","blue"),lty = c(2,1,3),lwd = 6,box.lty =0,cex =1.7)
fbarlab=substitute(bar(F)[X-Y],list(X=4,Y=7))
mtext(text=fbarlab,cex=1.5,side=2,line=-1,outer=TRUE)



par(oma=c(1,2,1,1))
catchplot(fit2,main = "Aggregated catch haddock model 2",cex.main = 2,cex.lab = 1.5,cex.axis=1.5,ylabInternal = "",ylim = c(0,400000))
mtext(text="Catch (tons)",cex=1.5,side=2,line=0.4,outer=TRUE)

par(oma=c(1,2,1,1))
catchplot(fit,main = "Aggregated catch haddock model 1",cex.main = 2,cex.lab = 1.5,cex.axis=1.5,ylabInternal = "",ylim = c(0,400000))
mtext(text="Catch (tons)",cex=1.5,side=2,line=0.4,outer=TRUE)

par(oma=c(1,2,1,1))
catchplot(fit3,main = "Aggregated catch haddock model 3",cex.main = 2,cex.lab = 1.5,cex.axis=1.5,ylabInternal = "",ylim = c(0,400000))
mtext(text="Catch (tons)",cex=1.5,side=2,line=0.4,outer=TRUE)





plot(res2,zmax = 2.5) #Needs small change in SAM to use zmax.
mtext("OSA residuals model 2 haddock", line = 1.5,at = 0, cex = 2.5)
plot(res,zmax = 2.5)
mtext("OSA residuals model 1 haddock", line = 1.5,at = 0, cex = 2.5)
plot(res3,zmax = 2.5)
mtext("OSA residuals model 3 haddock", line = 1.5,at = 0, cex = 2.5)


#Trompetplot
hh = which(res$fleet ==1 & res$age>4)
plot(res$mean[hh] ,res$residual[hh],ylim = c(-4,4.5),xlim = c(1,12),
     ylab = "OSA residual", xlab = "Predicted log-observation",
     main = "Catch of age 5+ haddock, model 1",
     cex.main = 1.8,cex.lab = 1.5,cex = 1.4,cex.axis = 1.4)
abline(h = 0)

hh = which(res2$fleet ==1 & res2$age>4)
plot(res2$mean[hh] ,res2$residual[hh],ylim = c(-4,4.5),xlim = c(1,12),
     ylab = "OSA residual", xlab = "Predicted log-observation",
     main = "Catch of age 5+ haddock, model 2",
     cex.main = 1.8,cex.lab = 1.5,cex = 1.4,cex.axis = 1.4)
abline(h = 0)



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

plotLink(fitSandard = fit,fitLink = fit2,f = 1,a = 1,main = "Catch sd")
plotLink(fitSandard = fit,fitLink = fit2,f = 1,a = 2,main = "Catch sd")
plotLink(fitSandard = fit,fitLink = fit2,f = 1,a = 3,main = "Catch sd")
plotLink(fitSandard = fit,fitLink = fit2,f = 2,a = 1,main = "Survey 1 sd")
plotLink(fitSandard = fit,fitLink = fit2,f = 3,a = 1,main = "Survey 2 sd")
plotLink(fitSandard = fit,fitLink = fit2,f = 4,a = 1,main = "Survey 3 sd")
plotLink(fitSandard = fit,fitLink = fit2,f = 5,a = 1,main = "Survey 4 sd")



