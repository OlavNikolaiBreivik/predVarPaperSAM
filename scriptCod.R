devtools::install_github("fishfollower/SAM/stockassessment",ref = "devPredVar")
library(stockassessment)

load("code/data/NEAcod/model.RData")
dat = fit$data
conf = fit$conf

conf<-loadConf(dat,"code/conf/NEAcod/modelLink.cfg", patch=TRUE)
par<-defpar(dat,conf)
fit2<-sam.fit(dat,conf,par)

conf<-loadConf(dat,"code/conf/NEAcod/modelFree.cfg", patch=TRUE)
par<-defpar(dat,conf)
fit3<-sam.fit(dat,conf,par)

AIC(fit,fit2,fit3)

#Parameter estimates
pl = as.list(fit2$sdrep,"Est")
plSd =  as.list(fit2$sdrep,"Std")
exp(pl$logSdLogObs)
exp(pl$logSdLogObs - 2*plSd$logSdLogObs)
exp(pl$logSdLogObs + 2*plSd$logSdLogObs)
exp(pl$predVarObs)
exp(pl$predVarObs - 2*plSd$predVarObs)
exp(pl$predVarObs + 2*plSd$predVarObs)


#jit
set.seed(123)
jj = jit(fit2,nojit = 20)

#Simulation
set.seed(123)
nsim = 20
sim = simstudy(fit2,nsim = nsim)

#OSA
res = residuals(fit)
res2 = residuals(fit2)
res3 = residuals(fit3)

#save(list = c("fit","fit2","fit3","jj","res","res2","res3","sim"),file = "resultsCod.RData")

nsim = length(sim)
par(oma=c(1,2,1,1))
ssbplot(fit2,main ="Simulation study cod",cex.main = 2.5,cex.lab = 1.5,cex.axis=1.3,ylabInternal = "",ci = FALSE)
for(i in 1:nsim){
  ssbplot(sim[[i]],ci = FALSE, add = TRUE, col = 'green')
}
ssbplot(fit2, ci = FALSE,add = TRUE)


#Boxplot
trueBeta = exp(fit2$pl$predVarObs)
trueAlpha = exp(fit2$pl$logSdLogObs)
estBeta = matrix(NA,nsim,length(exp(fit2$pl$predVarObs)))
estAlpha = matrix(NA,nsim,length(exp(fit2$pl$predVarObs)))
for(i in 1:nsim){
  estBeta[i,] = exp(sim[[i]]$pl$predVarObs)
  estAlpha[i,] = exp(sim[[i]]$pl$logSdLogObs)
}
boxplot(estBeta,xlab = "Beta-parameter",ylab = "Estimate",
        main = "Simulation experiment, cod, beta",cex.main = 2,
        cex.lab = 1.5,
        xaxt = 'n',cex.axis = 1.4)
names = c(expression(beta^("c")),
          expression(beta^("s1")),expression(beta^("s2")),expression(beta^("s3")),expression(beta^("s4")))
axis(side = 1, at=1:length(names), labels=names,cex.axis = 1.4)
points(trueBeta,col = 'red',cex = 2,pty = 3)

boxplot(estAlpha,xlab = "Alpha-parameter",ylab = "Estimate",
        main = "Simulation experiment, cod, alpha",ylim = c(0,80),cex.main =2,
        cex.lab = 1.5,
        xaxt = 'n',cex.axis = 1.4)
names = c(expression(alpha^("c")),expression(alpha[4]^("c")),expression(alpha[5-13]^("c")),
          expression(alpha^("s1")),expression(alpha^("s2")),expression(alpha^("s3")),expression(alpha^("s4")))
axis(side = 1, at=1:length(names), labels=names,cex.axis = 1.4)
points(trueAlpha,col = 'red',cex = 2,pty = 3)





par(oma=c(2,2,1,1))
ssbplot(fit2,main ="Cod SSB",cex.main = 2.5,cex.lab = 1.5,cex.axis=1.3)
ssb1 = exp(as.list(fit$sdrep,"Est",report = TRUE)$logssb)
lines(rev(fit$data$years), rev(ssb1), lwd=7, lty=2, col="red")
ssb3 = exp(as.list(fit3$sdrep,"Est",report = TRUE)$logssb)
lines(rev(fit$data$years), rev(ssb3), lwd=7, lty=3, col="blue")



par(oma=c(1,1.5,1,1))
fbarplot(fit2,partial = FALSE,main ="Cod fishing mortality",cex.main = 2.5,cex.lab = 1.3,cex.axis=1.5)
fbar1 = exp(as.list(fit$sdrep,"Est",report = TRUE)$logfbar)
fbar1 = fbar1[-length(fbar1)]
years = fit$data$years
years = years[-length(years)]
lines(rev(years), rev(fbar1), lwd=7, lty=2, col="red")
fbar3 = exp(as.list(fit3$sdrep,"Est",report = TRUE)$logfbar)
fbar3 = fbar3[-length(fbar3)]
lines(rev(years), rev(fbar3), lwd=7, lty=3, col="blue")


catchplot(fit,main = "Aggregated catch cod model 1",cex.main = 2,cex.lab = 1.5,cex.axis=1.5)

catchplot(fit2,main = "Aggregated catch cod model 2",cex.main = 2,cex.lab = 1.5,cex.axis=1.5)

catchplot(fit3,main = "Aggregated catch cod model 3",cex.main = 2,cex.lab = 1.5,cex.axis=1.5)



#Plot only catch residuals
plotRes = function(x,zmax,...){
  class(x)<-NULL
  x = as.data.frame(x)
  x = x[x$fleet==1,]
  class(x)<-"samres"

  add_legend <- function(x,zmax = NULL,...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n',...)
    zscale <- pretty(x$residual,min.n=4)
    uu<-par("usr")
    yy<-rep(uu[3]+.03*(uu[4]-uu[3]), length(zscale))
    xx<-seq(uu[1]+.10*(uu[2]-uu[1]),uu[1]+.4*(uu[2]-uu[1]), length=length(zscale))
    text(xx,yy,labels=zscale)
    colb <- ifelse(zscale<0, rgb(1, 0, 0, alpha=.5), rgb(0, 0, 1, alpha=.5))
    bs<-1
    if(is.null(zmax)){
      points(xx,yy,cex=sqrt(abs(zscale))/max(sqrt(abs(zscale)), na.rm=TRUE)*5*bs, pch=19, col=colb)
    }else{
      points(xx,yy,cex=sqrt(abs(zscale))/zmax*5*bs, pch=19, col=colb)
    }
  }
  plotby(x$year, x$age, x$residual, by="Residual catch", xlab="", ylab="",zmax = zmax,...)
  add_legend(x, zmax = zmax)
}


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

plotLink = function(fitSandard,fitLink,f,aa,col = 1,main = "",ylim = c(0,3)){

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
ylim = c(0,2.5)
plotLink(fitSandard = fit,fitLink = fit2,f = 1,aa = 1,main = "Catch sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 2,aa = 1,main = "Survey 1 sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 3,aa = 1,main = "Survey 2 sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 4,aa = 1,main = "Survey 3 sd",ylim = ylim)
plotLink(fitSandard = fit,fitLink = fit2,f = 5,aa = 1,main = "Survey 4 sd",ylim = ylim)


#Assign catchability from model 1 to model 2
conf<-loadConf(dat,"code/conf/NEAcod/modelLink.cfg", patch=TRUE)
par<-defpar(dat,conf)
par$logFpar = fit$pl$logFpar
map = list()
map$logFpar = as.factor(rep(NA,length(fit$pl$logFpar)))
fitTest<-sam.fit(dat,conf,par,map = map)

ssbplot(c(fitTest,fit),main ="SSB with fixed catchability",cex.main = 2.5,cex.lab = 1.5,cex.axis=1.3)
