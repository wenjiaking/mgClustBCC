#' Log likelihood by plugging in posterior samples
#'
#' @param X N by G omics matrix, where each row is a sample
#' @param Y N by O outcome matrix, where each row is a sample
#' @param Z covariates in the outcome association model, N by q matrix without column of intercept
#' @param MCMCres output of MCMC function
#' @param nsamples the number of posterior samples after burn-in (must smaller than n_iter in MCMC function), 1000 by default
#' @param Y.ind censor indicators. NULL if none of the outcomes in Y is survival, otherwise N by O matrix where non-survival outcomes are NA
#'
#' @return
#' the median of log-likelihood over all posterior samples
#' @export
#'
#' @examples
#' \dontrun{
#' test.loglkh=postLogLK(X=LungData$X[,1:1000],Y=LungData$Y[,1:3],
#' Y.ind=NULL,Z=LungData$Z[,-1],nsamples=500,MCMCres=test)
#' }
postLogLK=function(X,Y,Z,MCMCres,nsamples=1000,Y.ind=NULL) {
  ntotal=length(MCMCres$BETA0.steps)
  nstart=ntotal-nsamples+1
  loglkh_vec=sapply(nstart:ntotal, function(i) {
    posterior_s=list(BETA0.est=MCMCres$BETA0.steps[[i]],
                     BETA.est=MCMCres$BETA.steps[[i]],
                     sigma.est=MCMCres$pars.steps[i,grep("sigma",colnames(MCMCres$pars.steps)),drop=T],
                     s.est=MCMCres$s.steps[i,],
                     TAU.est=MCMCres$pars.steps[i,grep("tau",colnames(MCMCres$pars.steps))],
                     MU.est=MCMCres$MU.steps[[i]],
                     C.est=MCMCres$C.steps[i,],
                     L_G.est=MCMCres$L_G.steps[i,],
                     L_O.est=MCMCres$L_O.steps[[i]],
                     alpha.est=MCMCres$pars.steps[i,grep("alpha",colnames(MCMCres$pars.steps))])
    if (sum(posterior_s$s.est)<20) {
      loglkh.tmp=NA
    }
    else {
      loglkh.tmp=loglkh(X=X,Y=Y,Z=Z,estimates=posterior_s,Y.ind=Y.ind)
    }
    return(loglkh.tmp)
  })
  return(median(loglkh_vec,na.rm = T))
}
loglkh=function(X,Y,Z,estimates,Y.ind=NULL) {
  BETA0.est=estimates$BETA0.est
  BETA.est=estimates$BETA.est
  sigma.est=estimates$sigma.est
  s.est=estimates$s.est
  TAU.est=estimates$TAU.est[(s.est==1)] #if the different chains selected very different number of genes, the omics likelihood may be in different scale!!
  MU.est=estimates$MU.est[(s.est==1),] #but cannot calculate the density of full-dimension of genes, which is too high dimension resulting errors
  C.est=as.numeric(estimates$C.est)
  L_G.est=as.numeric(estimates$L_G.est)
  L_O.est=estimates$L_O.est
  alpha.est=estimates$alpha.est
  X=scale(X,center = T,scale = F)
  X=X[,(s.est==1)]
  lkh.omics=sum(sapply(1:nrow(X),function(i) {log(dmvnorm(X[i,,drop=F],MU.est[,L_G.est[i]],diag(TAU.est)))}))
  # lkh.omics=sum(sapply(1:ncol(X),function(g) log(dnorm(X[i,g],mean=MU.est[g,L_G.est[i]],sd=sqrt(TAU.est[g])))))

  n_O=ncol(Y)
  for (o in 1:n_O) {
    BETA0_outcomes=BETA0.est[o,as.numeric(L_O.est[,o])]
    BETA=BETA.est[o,]
    y=as.vector(Y[,o])
    sigma=sigma.est[o]
    if (!is.null(Y.ind)) {
      y.ind=as.vector(Y.ind[,o])
    }
    else {
      y.ind=rep(NA,length(y))
    }

    if (all(is.na(y.ind))) {
      if (length(table(y))==2) {
        lkh.outcome=sum((BETA0_outcomes+Z%*%BETA)*y-log(1+exp(BETA0_outcomes+Z%*%BETA)))
      }
      else {
        lkh.outcome=sum(sapply(1:length(y),function(i) {log(dnorm(y[i],BETA0_outcomes[i]+sum(BETA*Z[i,]),sigma))}))
      }
    }
    else {
      lkh.outcome=sum(1/sqrt(sigma)*(y-BETA0_outcomes-Z%*%BETA)-log(1+exp(1/sqrt(sigma)*(y-BETA0_outcomes-Z%*%BETA)))*(1+y.ind))
    }
  }

  # if (ncol(Y)>1) {
  #   BETA0_outcomes=sapply(1:ncol(L_O.est),function(i) BETA0.est[i,as.numeric(L_O.est[,i])])
  #   lkh.outcome=sum(sapply(1:nrow(Y),function(i) {log(dmvnorm(Y[i,,drop=F],BETA0_outcomes[i,]+BETA.est%*%Z[i,],diag(sigma.est)))}))
  # }
  # else {
  #   lkh.outcome=sum(log(dnorm(as.numeric(Y),as.numeric(BETA0.est[as.numeric(L_O.est)]),sigma.est)))
  # }

  K=ncol(MU.est)
  lkh.L=rep(alpha.est[1],nrow(X))
  lkh.L[C.est!=L_G.est]=(1-alpha.est[1])/(K-1)
  lkh.L=sum(log(lkh.L))

  # nO=sum(alpha.est[-1]>alpha_O.cut)
  n_O=ncol(Y)
  lkh.O=sum(sapply(1:n_O,function(o) {
    lkh.O.tmp=rep(alpha.est[1+o],nrow(X))
    lkh.O.tmp[C.est!=as.numeric(L_O.est[,o])]=(1-alpha.est[1+o])/(K-1)
    lkh.O.tmp=sum(log(lkh.O.tmp))
  }))

  return(lkh.omics+lkh.outcome+lkh.L+lkh.O)
}
