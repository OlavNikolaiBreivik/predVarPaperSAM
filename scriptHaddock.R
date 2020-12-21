devtools::install_github("fishfollower/SAM/stockassessment",ref = "devPredVar")
library(stockassessment)

#Currently used fit in assessment
load("code/data/NEAhaddock/model.RData")
dat = fit$data
conf = fit$conf
par<-defpar(dat,conf)

fit <- sam.fit(dat,conf,par)

conf = fit$conf
conf$predVarObsLink = matrix(-1, 5, 11)
par<-defpar(dat,conf)
fit2 <- sam.fit(dat,conf,par)

conf3 = loadConf(dat,file = "code/conf/NEAhaddock/modelFree.cfg")
par3<-defpar(dat,conf3)
fit3 <- sam.fit(dat,conf3,par3)


AIC(fit,fit2,fit3)

#jit
set.seed(123)
jj = jit(fit,nojit = 20)

#OSA
res = residuals(fit)
res2 = residuals(fit2)
res3 = residuals(fit3)

#sim
set.seed(123)
nsim = 20
sim = simstudy(fit,nsim = nsim)

#save(list = c("fit","fit2","fit3","jj","res","res2","res3","sim"),file = "resultsHad.RData")


#Parameter estimates
pl = as.list(fit$sdrep,"Est")
plSd =  as.list(fit$sdrep,"Std")
exp(pl$logSdLogObs)
exp(pl$logSdLogObs - 2*plSd$logSdLogObs)
exp(pl$logSdLogObs + 2*plSd$logSdLogObs)
exp(pl$predVarObs)
exp(pl$predVarObs - 2*plSd$predVarObs)
exp(pl$predVarObs + 2*plSd$predVarObs)


ssbplot(fit,main ="Haddock SSB",cex.main = 2.5,cex.lab = 1.5,cex.axis=1.3)
ssb1 = exp(as.list(fit2$sdrep,"Est",report = TRUE)$logssb)
lines(rev(fit$data$years), rev(ssb1), lwd=7, lty=2, col="red")
ssb3 = exp(as.list(fit3$sdrep,"Est",report = TRUE)$logssb)
lines(rev(fit$data$years), rev(ssb3), lwd=7, lty=3, col="blue")
legend(x=1950,y=650000,legend = c("Model 1", "Model 2", "Model 3")
       ,col=c("red","black","blue"),lty = c(2,1,3),lwd = 6,box.lty =0,cex =1.7)



fbarplot(fit,partial = FALSE,main ="Haddock fishing mortality",cex.main = 2.5,cex.lab = 1.3,cex.axis=1.5)
fbar1 = exp(as.list(fit2$sdrep,"Est",report = TRUE)$logfbar)
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


#Catch plot
catchplot(fit,main = "Aggregated catch haddock model 2",cex.main = 2,cex.lab = 1.5,cex.axis=1.5)
catchplot(fit2,main = "Aggregated catch haddock model 1",cex.main = 2,cex.lab = 1.5,cex.axis=1.5)
catchplot(fit3,main = "Aggregated catch haddock model 3",cex.main = 2,cex.lab = 1.5,cex.axis=1.5)





plot(res,zmax = 2.5) #Needs small change in SAM to use zmax.
mtext("OSA residuals model 2 haddock", line = 1.5,at = 0, cex = 2.5)
plot(res2,zmax = 2.5)
mtext("OSA residuals model 1 haddock", line = 1.5,at = 0, cex = 2.5)
plot(res3,zmax = 2.5)
mtext("OSA residuals model 3 haddock", line = 1.5,at = 0, cex = 2.5)


#Trompetplot
hh = which(res2$fleet ==1 & res2$age>4)
plot(res2$mean[hh] ,res2$residual[hh],ylim = c(-4,4.5),xlim = c(1,12),
     ylab = "OSA residual", xlab = "Predicted log-observation",
     main = "Catch of age 5+ haddock, model 1",
     cex.main = 1.8,cex.lab = 1.5,cex = 1.4,cex.axis = 1.4)
abline(h = 0)

hh = which(res$fleet ==1 & res$age>4)
plot(res$mean[hh] ,res$residual[hh],ylim = c(-4,4.5),xlim = c(1,12),
     ylab = "OSA residual", xlab = "Predicted log-observation",
     main = "Catch of age 5+ haddock, model 2",
     cex.main = 1.8,cex.lab = 1.5,cex = 1.4,cex.axis = 1.4)
abline(h = 0)





ssbplot(fit,main ="Haddock SSB",cex.main = 2.5,cex.lab = 1.5,cex.axis=1.3,ci = FALSE)
for(i in 1:nsim){
  ssbplot(sim[[i]],ci = FALSE, add = TRUE, col = 'green')
}
ssbplot(fit, ci = FALSE,add = TRUE)

#Boxplot
trueBeta = exp(fit$pl$predVarObs)
trueAlpha = exp(fit$pl$logSdLogObs)
estBeta = matrix(NA,nsim,length(exp(fit$pl$predVarObs)))
estAlpha = matrix(NA,nsim,length(exp(fit$pl$predVarObs)))
for(i in 1:nsim){
  estBeta[i,] = exp(sim[[i]]$pl$predVarObs)
  estAlpha[i,] = exp(sim[[i]]$pl$logSdLogObs)
}
boxplot(estBeta,xlab = "Beta-parameter",ylab = "Estimate",
        main = "Simulation experiment, haddock, beta",cex.main = 2,
        cex.lab = 1.5,
        xaxt = 'n',cex.axis = 1.4)
names = c(expression(beta^("c")),expression(beta[4]^("c")),expression(beta[5-13]^("c")),
          expression(beta^("s1")),expression(beta^("s2")),expression(beta^("s3")),expression(beta^("s4")))
axis(side = 1, at=1:length(names), labels=names,cex.axis = 1.4)
points(trueBeta,col = 'red',cex = 2,pty = 3)

boxplot(estAlpha,xlab = "Alpha-parameter",ylab = "Estimate",
        main = "Simulation experiment, haddock, alpha",ylim = c(0,10),cex.main = 2,
        cex.lab = 1.5,
        xaxt = 'n',cex.axis = 1.4)
names = c(expression(alpha^("c")),expression(alpha[4]^("c")),expression(alpha[5-13]^("c")),
          expression(alpha^("s1")),expression(alpha^("s2")),expression(alpha^("s3")),expression(alpha^("s4")))
axis(side = 1, at=1:length(names), labels=names,cex.axis = 1.4)
points(trueAlpha,col = 'red',cex = 2,pty = 3)





#Plot of sd
N = 20000
ss = rmvnorm(N,mu = fit$opt$par,Sigma = fit$sdrep$cov.fixed)
alpha = exp(ss[,which(names(fit$opt$par)=="logSdLogObs")])
beta = exp(ss[,which(names(fit$opt$par)=="predVarObs")])

plotLink = function(fitSandard,fitLink,f,aa,col = 1,main = "",ylim =c(0,3)){

  if(f>1){
    aa = aa + max(fitLink$conf$keyVarObs[f-1,])+1
  }
  age = which(fitLink$conf$keyVarObs[f,]==aa-1) + fitLink$conf$minAge -1
  obsIndex = which(fitLink$data$aux[,2]==f & fitLink$data$aux[,3] %in% age)
  obs = fitLink$data$logobs[obsIndex]
  obs = obs[!is.na(obs)]
  l = seq(min(obs),max(obs), by = 0.01)

  sd = matrix(0,N,length(l))
  for(i in 1:N){
    sd[i,] = sqrt(log(alpha[i,aa]* exp(l*(beta[i,aa]-1)) +1))
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
         ylab = "Sd", xlab = "Log-prediction",cex.main = 1.7,
         cex.lab = 1.5,
         cex.axis = 1.2)
  }else{
    plot(l,sdM,type = 'l', ylim = ylim,lwd = 2,col = col,main = paste0(main, " age ",min(age),"-",max(age)),
         ylab = "Sd", xlab = "Log-prediction",cex.main = 1.7,
         cex.lab = 1.5,
         cex.axis = 1.2)
  }
  lines(l,sdL,col = col)
  lines(l,sdU,col = col)
  abline(h = exp(fitSandard$pl$logSdLogObs[aa]), col = col)
  abline(h = exp(fitSandard$pl$logSdLogObs[aa] + 1.96*fitSandard$plsd$logSdLogObs[aa]),col = col)
  abline(h = exp(fitSandard$pl$logSdLogObs[aa] - 1.96*fitSandard$plsd$logSdLogObs[aa]),col = col)

  points(obs,rep(exp(fitSandard$pl$logSdLogObs[aa]),length(obs)))
}
ylim = c(0,1.7)
plotLink(fitSandard = fit2,fitLink = fit,f = 1,aa = 1,main = "Catch sd")
plotLink(fitSandard = fit2,fitLink = fit,f = 1,aa = 2,main = "Catch sd")
plotLink(fitSandard = fit2,fitLink = fit,f = 1,aa = 3,main = "Catch sd")
plotLink(fitSandard = fit2,fitLink = fit,f = 2,aa = 1,main = "Survey 1 sd")
plotLink(fitSandard = fit2,fitLink = fit,f = 3,aa = 1,main = "Survey 2 sd")
plotLink(fitSandard = fit2,fitLink = fit,f = 4,aa = 1,main = "Survey 3 sd")
plotLink(fitSandard = fit2,fitLink = fit,f = 5,aa = 1,main = "Survey 4 sd")



