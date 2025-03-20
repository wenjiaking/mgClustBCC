#' Predict the cluster labels
#'
#' @param X N by G omics matrix, where each row is a sample
#' @param Y N by O outcome matrix, where each row is a sample (optional)
#' @param Z covariates in the outcome association model, N by q matrix without column of intercept
#' @param estimates output of estimate function
#' @param K the number of clusters
#' @param Y.ind censor indicators. NULL if none of the outcomes in Y is survival, otherwise N by O matrix where non-survival outcomes are NA
#'
#' @return
#' A list of predicted cluster labels.
#' \itemize{
#' \item{C.pred: }{a vector of length N, the predicted master cluster label for each sample}
#' \item{L_G.pred: }{a vector of length N, the predicted omics-specific cluster label for each sample}
#' \item{L_O.pred: }{a N by O matrix where each row represents the predicted outcome-specific cluster label for each sample}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' test.pred=MCMCpredict(X=LungData$X[,1:1000],Y=LungData$Y[,1:3],
#' Y.ind=NULL,Z=LungData$Z[,-1],K=3,estimates=test.est)
#' }
MCMCpredict=function(X=X,Y=Y,Z=Z,estimates,K,Y.ind=NULL) {
  X=scale(X,center = T,scale = F)
  X=X[,(estimates$s.est==1)]
  MU.est=estimates$MU.est[(estimates$s.est==1),]
  TAU.est=estimates$TAU.est[(estimates$s.est==1)]
  D.test=t(X)
  mult_pdf<-matrix(,nrow=ncol(D.test),ncol=ncol(MU.est))
  for(j1 in 1:ncol(MU.est)){
    temp<-mult_density1(as.numeric(D.test),mu=rep(MU.est[,j1],times=ncol(D.test)),sigma=rep(TAU.est,times=ncol(D.test)))
    temp_matrix<-matrix(temp,nrow=nrow(MU.est),ncol=ncol(D.test),byrow=FALSE)
    mult_pdf[,j1]<-t(temp_matrix)%*%rep(1,nrow(MU.est))
  }
  L_G.pred=apply(mult_pdf,1,which.max)
  if (!is.null(Y)) {
    BETA.est=estimates$BETA.est
    BETA0.est=estimates$BETA0.est
    sigma.est=estimates$sigma.est
    n_O=ncol(Y)
    N=nrow(Y)
    L_O.pred=matrix(,,nrow=N,ncol=n_O)
    for (o in 1:n_O) {
      y=Y[,o]
      beta0.est=matrix(rep(BETA0.est[o,],each=N),nrow=N)
      if (!is.null(Y.ind)) {
        y.ind=as.vector(Y.ind[,o])
      }
      else {
        y.ind=rep(NA,length(y))
      }

      if (all(is.na(y.ind))) {
        if (length(table(y))==2) {
          y.density=-(beta0.est+as.vector(Z %*% BETA.est[o,]))*y+log(1+exp(beta0.est+as.vector(Z %*% BETA.est[o,])))
        }
        else {
          y.density=(beta0.est+as.vector(Z %*% BETA.est[o,]-y))^2 #loss
        }
      }
      else {
        sigma=sigma.est[o]
        y.density=-1/sqrt(sigma)*(beta0.est+as.vector(Z %*% BETA.est[o,]-y))+log(1+exp(1/sqrt(sigma)*(beta0.est+as.vector(Z %*% BETA.est[o,]-y))))*(1+y.ind)
      }

      O.pred=apply(y.density,1,which.min)
      L_O.pred[,o]=O.pred
    }

    pi.est=matrix(rep(estimates$pi.est,each=N),N)
    alpha.est=estimates$alpha.est
    alphaG.est=alpha.est["alpha_G"]
    probG=t(matrix(rep((1-alphaG.est)/(K-1),K*N),nrow=N))
    probG[(L_G.pred+K*(0:(N-1)))]=alphaG.est
    probC=pi.est*t(probG)
    for (o in 1:n_O) {
      alphaO.est=alpha.est[paste0("alpha_O_",o)]
      probO=t(matrix(rep((1-alphaO.est)/(K-1),K*N),nrow=N))
      O.pred=L_O.pred[,o]
      probO[(O.pred+K*(0:(N-1)))]=alphaO.est
      probC=probC*t(probO)
    }
    C.pred=apply(probC,1,which.max)
  }
  else {
    L_O.pred=NULL
    C.pred=NULL
  }

  return(list(C.pred=C.pred,L_G.pred=L_G.pred,L_O.pred=L_O.pred))
}

mult_density1<-function(x,mu,sigma){
  lik<-dnorm(x,mean=mu,sd=sqrt(sigma),log=TRUE)
  return(lik)
}
