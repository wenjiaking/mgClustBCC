#' Point estimates based on posterior samples
#'
#' @param MCMCres output of MCMC function
#' @param nsamples the number of posterior samples after burn-in (must smaller than n_iter in MCMC function), 1000 by default
#' @param K the number of clusters
#'
#' @return
#' A list of point estimates.
#' \itemize{
#' \item{s.est: }{a vector of length G, indicating whether the feature is selected (1) or not (0)}
#' \item{alpha.est: }{a vector contains the estiamtes of omics-specific consensus rate and outcome-specific consensus rates}
#' \item{pi.est: }{a vector of length K, the estiamtes of cluster probabilities}
#' \item{MU.est: }{a G by K matrix where each row represents the omics-specific cluster centers for each feature}
#' \item{BETA0.est: }{an O by K matrix where each row represents the intercepts across outcome-specific cluster for each outcome}
#' \item{BETA.est: }{an O by q matrix where each row represents the covariate coefficients for each outcome}
#' \item{C.est: }{a vector of length N, the master cluster estimate for each sample}
#' \item{L_G.est: }{a vector of length N, the omics-specific cluster estimate for each sample}
#' \item{L_O.est: }{a N by O matrix where each row represents the outcome-specific cluster estimates for each sample}
#' \item{TAU.est: }{a vector of length G, the estimates of omics variance}
#' \item{sigma.est: }{a vector of length O, the estimates of outcome variance}
#' \item{e.est: }{a vector of length 2, the variance estiamtes of cluster-varying and stable feature}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   test.est=estimate(test,nsamples = 500,K=3)
#'}
estimate=function(MCMCres,nsamples=1000,K) {
  s.steps=MCMCres$s.steps
  ntotal=nrow(s.steps)
  nstart=ntotal-nsamples+1
  s.steps=s.steps[nstart:ntotal,]
  s.mod=apply(s.steps,2,function(x) names(which.max(table(x))))
  MU.steps=MCMCres$MU.steps[nstart:ntotal]
  G=nrow(MU.steps[[1]])
  # index=which(s.mod==1)
  # MU.est=matrix(,,nrow=length(index),ncol=K)
  # for (k in 1:K) {
  #   MU.est[,k]=sapply(index+G*(k-1),function(i) median(sapply(MU.steps,"[[",i)))
  # }

  MU.est=matrix(,,nrow=G,ncol=K)
  for (k in 1:K) {
    MU.est[,k]=sapply((1:G)+G*(k-1),function(i) median(sapply(MU.steps,"[[",i)))
  }
  TAU.steps=MCMCres$pars.steps[nstart:ntotal,grep("tau",colnames(MCMCres$pars.steps))]
  # TAU.est=apply(TAU.steps[,s.mod==1],2,median)
  TAU.est=apply(TAU.steps,2,median)
  BETA0.steps=MCMCres$BETA0.steps[nstart:ntotal]
  n_O=nrow(BETA0.steps[[1]])
  BETA0.est=matrix(sapply(1:(K*n_O),function(i) median(sapply(BETA0.steps,"[[",i))),nrow=n_O)
  BETA.steps=MCMCres$BETA.steps[nstart:ntotal]
  q=ncol(BETA.steps[[1]])
  BETA.est=matrix(sapply(1:(q*n_O),function(i) median(sapply(BETA.steps,"[[",i))),nrow=n_O)

  C.steps=MCMCres$C.steps[nstart:ntotal,]
  C.mod=apply(C.steps,2,function(x) names(which.max(table(x))))
  L_G.steps=MCMCres$L_G.steps[nstart:ntotal,]
  L_G.mode=apply(L_G.steps,2,function(x) names(which.max(table(x))))

  L_O.steps=MCMCres$L_O.steps[nstart:ntotal]
  N=nrow(L_O.steps[[1]])
  L_O.mod=matrix(sapply(1:(n_O*N), function(i) {
    names(which.max(table(sapply(L_O.steps,"[[",i))))
  }),ncol = n_O)

  parsmat=MCMCres$pars.steps[nstart:ntotal,]
  alphamat=parsmat[,grep("alpha",colnames(parsmat))]
  alpha.est=apply(alphamat,2,median)
  pimat=parsmat[,grep("pi",colnames(parsmat))]
  pi.est=apply(pimat,2,median)

  sigma.steps=parsmat[,grep("sigma",colnames(parsmat)),drop=F]
  sigma.est=apply(sigma.steps,2,median)

  e.steps=parsmat[,c("e1","e0")]
  e.est=apply(e.steps,2,median)

  return(list(s.est=s.mod,alpha.est=alpha.est,pi.est=pi.est,MU.est=MU.est,BETA0.est=BETA0.est,BETA.est=BETA.est,
              C.est=C.mod,L_G.est=L_G.mode,L_O.est=L_O.mod,TAU.est=TAU.est,sigma.est=sigma.est,e.est=e.est))
}
