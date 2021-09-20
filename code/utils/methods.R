
#sim
simStudyPaper = function(simdata,fit, nsim, ncores = detectCores(), newtonsteps=0,...)
{
  if (ncores > 1) {
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl))
    lib.ver <- dirname(path.package("stockassessment"))
    clusterExport(cl, varlist = c("simdata","fit", "lib.ver", "newtonsteps"), envir = environment())
    clusterEvalQ(cl, {
      library(stockassessment, lib.loc = lib.ver)
    })
    runs <- parLapply(cl, 1:nsim, function(i) sam.fit(simdata[[i]],
                                                      fit$conf, fit$obj$env$par, newtonsteps = newtonsteps))
  }
  else {
    runs <- lapply(1:nsim, function(i) sam.fit(simdata[[i]], fit$conf,
                                               fit$obj$env$par, newtonsteps = newtonsteps))
  }
  class(runs) <- "samset"
  runs
}
