#' Metropolis-within-Gibbs Sampler for mgClustBCC Inference
#'
#' @param n_iter the number of iteration
#' @param K the number of clusters
#' @param X N by G omics matrix, where each row is a sample
#' @param Y N by O outcome matrix, where each row is a sample
#' @param Y.ind censor indicators. NULL if none of the outcomes in Y is survival, otherwise N by O matrix where non-survival outcomes are NA
#' @param Z covariates in the outcome association model, N by (q+1) matrix where the first column is the intercept of value 1.
#' @param alpha_G.init initial value of omics-specific consensus rate, 0.5 by default
#' @param alpha_O.init initial value of outcome-specific consensus rates, a vector of the same length as the number of outcomes, 0.5 by default
#' @param v0_vec variance scale of cluster-specific intercepts, a vector of the same length as the number of outcomes, 100 by default
#' @param delta step size of Metropolis-Hasting sampler for survival outcome, 0.5 by default
#' @param d1 shape of inverse gamma prior for omics variance tau, 0.001 by default
#' @param d2 rate of inverse gamma prior for omics variance tau, 0.001 by default
#' @param t1 shape of inverse gamma prior for outcome variance sigma (continous or survival), 2 by default
#' @param t2 rate of inverse gamma prior for outcome variance sigma (continous or survival), 2 by default
#' @param a1 shape of inverse gamma prior for the variance of cluster-varying genes, 4 by default
#' @param a2 rate of inverse gamma prior for the variance of cluster-varying genes, 450 by default
#' @param b1 shape of inverse gamma prior for the variance of stable genes, 2 by default
#' @param b2 rate of inverse gamma prior for the variance of stable genes, 0.005 by default
#' @param r_G1 shape1 of truncated beta prior for the omics-specific consensus rate, 5 by default
#' @param r_G2 shape2 of truncated beta prior for the omics-specific consensus rate, 1 by default
#' @param r_O1 shape1 of truncated beta prior for the outcome-specific consensus rates, a vector of the same length as the number of outcomes, 1 by default
#' @param r_O2 shape2 of truncated beta prior for the outcome-specific consensus rates, a vector of the same length as the number of outcomes, 1 by default
#' @param h parameters of Dirichlet prior, rep(1,K) by default
#' @param cutoff_G truncation point of truncated beta prior for the omics-specific consensus rate, 1/K by default
#' @param cutoff_O truncation point of truncated beta prior for the outcome-specific consensus rates, a vector of the same length as the number of outcomes, 1/K by default
#' @param init.p.adjust pre-screening based on original p-values or adjusted p-values, TRUE by default
#' @param init.pcut pre-screening p-value cutoff, 0.05 by default
#' @param init.NG pre-specified number of initially selected genes, NULL by default
#'
#' @return
#' A list of posterior samples.
#' \itemize{
#' \item{s.steps: }{n_iter by G binary matrix where each row represents the gene selection values at one iteration}
#' \item{C.steps: }{n_iter by N matrix of master cluster where each row represents master cluster labels of all samples at one iteration}
#' \item{L_G.steps: }{n_iter by N matrix of omics-specific cluster where each row represents omics-specific cluster labels of all samples at one iteration}
#' \item{L_O.steps: }{a list of n_iter elements, where each element is an N by O matrix of omics-specific clusters with each row representing cluster labels of one sample for all outcomes}
#' \item{BETA0.steps: }{a list of n_iter elements, where each element is an O by K matrix of intercepts across outcome-specific cluster for each outcome}
#' \item{BETA.steps: }{a list of n_iter elements, where each element is an O by q matrix of covariate coefficients for each outcome}
#' \item{MU.steps: }{a list of n_iter elements, where each element is a G by K matrix with each row representing the omics-specific cluster centers for each feature}
#' \item{pars.steps: }{a matrix of parameters, where each row represents omics-specific consensus rate, outcome-specific consensus rates, variance of cluster-varying and stable genes, omics variance, outcome variance, cluster probability, probablity of gene selection for one iteration}
#' }
#'
#' @export
#' @import stats survival MASS mclust BayesLogit
#'
#' @examples
#' \dontrun{
#'   data('LungData') #load dataset
#'   test=MCMC(n_iter=1000,K=K,X=LungData$X[,1:1000],Y=LungData$Y[,1:3],
#'   Y.ind=NULL,Z=LungData$Z)
#'}

MCMC=function(n_iter,K,X,Y,Y.ind=NULL,Z,alpha_G.init=0.5,alpha_O.init=0.5,v0_vec=100,delta=0.5,
              d1=0.001,d2=0.001,t1=2,t2=2,a1=4,a2=450,b1=2,b2=0.005,r_G1=5,r_G2=1,r_O1=1,r_O2=1,
              h=NULL,cutoff_G=NULL,cutoff_O=NULL,init.p.adjust=TRUE,init.pcut=0.05,init.NG=NULL) {
  N=nrow(Y)
  n_O=ncol(Y)
  G=ncol(X)
  q=ncol(Z)
  #standardize
  X=scale(X,center = T,scale=F)
  if (is.null(h)) {h=rep(1,K)}
  if (is.null(cutoff_G)) {cutoff_G=1/K}
  if (is.null(cutoff_O)) {cutoff_O=1/K}
  if (length(r_O1)<n_O) {r_O1=rep(r_O1,length.out=n_O)}
  if (length(r_O2)<n_O) {r_O2=rep(r_O2,length.out=n_O)}
  if (length(v0_vec)<n_O) {v0_vec=rep(v0_vec,length.out=n_O)}
  if (length(cutoff_O)<n_O) {cutoff_O=rep(cutoff_O,length.out=n_O)}
  if (length(alpha_O.init)<n_O) {alpha_O.init=rep(alpha_O.init,length.out=n_O)}

  sigma_vec=c()
  beta0_vec=c() #Maybe center Y and set betaO_vec to be 0?
  #v0_vec=c()
  beta_mat=matrix(NA,nrow=n_O,ncol=q-1)
  v_mat=matrix(NA,nrow=n_O,ncol=q-1)

  if (is.null(Y.ind)) {
    for (o in 1:n_O) {
      y=as.vector(Y[,o])
      if (length(table(y))==2) {
        #binary outcome, fit the logistic regression
        fit.O=glm(y~0+.,data=as.data.frame(cbind(y,Z)),family = binomial(link = "logit"))
        sigma_vec=c(sigma_vec,var(exp(predict(fit.O))/(1+exp(predict(fit.O)))-y))
        beta_temp=fit.O$coefficients
        beta0_vec=c(beta0_vec,beta_temp[1])
        beta_mat[o,]=beta_temp[-1]
        v_temp=(summary(fit.O)$coefficients[,2]^2)/sigma_vec[o]
        v_mat[o,]=v_temp[-1]
      }
      else {
        #continuous outcome, fit the linear regression
        fit.O=lm(y~0+.,data=as.data.frame(cbind(y,Z)))
        sigma_vec=c(sigma_vec,var(predict(fit.O)-y))
        beta_temp=fit.O$coefficients
        beta0_vec=c(beta0_vec,beta_temp[1])
        beta_mat[o,]=beta_temp[-1]
        v_temp=(summary(fit.O)$coefficients[,2]^2)/sigma_vec[o] #some covariates may not be cluster-varying
        #v0_vec=c(v0_vec,v_temp[1])
        v_mat[o,]=v_temp[-1]
      }
    }
  }
  else {
    for (o in 1:n_O) {
      y=as.vector(Y[,o])
      y.ind=as.vector(Y.ind[,o])
      if (all(is.na(y.ind))) {
        if (length(table(y))==2) {
          #binary outcome, fit the logistic regression
          fit.O=glm(y~0+.,data=as.data.frame(cbind(y,Z)),family = binomial(link = "logit"))
          sigma_vec=c(sigma_vec,var(exp(predict(fit.O))/(1+exp(predict(fit.O)))-y))
          beta_temp=fit.O$coefficients
          beta0_vec=c(beta0_vec,beta_temp[1])
          beta_mat[o,]=beta_temp[-1]
          v_temp=(summary(fit.O)$coefficients[,2]^2)/sigma_vec[o]
          v_mat[o,]=v_temp[-1]
        }
        else {
          #continuous outcome, fit the linear regression
          fit.O=lm(y~0+.,data=as.data.frame(cbind(y,Z)))
          sigma_vec=c(sigma_vec,var(predict(fit.O)-y))
          beta_temp=fit.O$coefficients
          beta0_vec=c(beta0_vec,beta_temp[1])
          beta_mat[o,]=beta_temp[-1]
          v_temp=(summary(fit.O)$coefficients[,2]^2)/sigma_vec[o] #some covariates may not be cluster-varying
          #v0_vec=c(v0_vec,v_temp[1])
          v_mat[o,]=v_temp[-1]
        }
      } else {
        #survival outcome, fit AFT model with log-logistic parametric family
        fit.O<-survreg(Surv(y, y.ind) ~ 0+., as.data.frame(cbind(y,y.ind,Z)),robust=T,
                       dist="loglogistic")
        sigma_vec=c(sigma_vec,fit.O$scale^2)
        beta_temp=fit.O$coefficients
        beta0_vec=c(beta0_vec,beta_temp[1])
        beta_mat[o,]=beta_temp[-1]
        v_temp=diag((summary(fit.O)$var))/sigma_vec[o]
        v_mat[o,]=v_temp[-c(1,length(v_temp))]
      }
    }
  }


  t.stat<-matrix(NA,nrow=n_O,ncol=ncol(X))
  if (is.null(Y.ind)) {
    for (o in 1:n_O) {
      y=as.vector(Y[,o])
      if (length(table(y))==2) {
        #binary outcome, fit the logistic regression
        for(i in 1:ncol(X)){
          data<-data.frame(y=y,x=X[,i])
          mod<-glm(y~x,data=data,family = binomial(link = "logit"))
          mod<-summary(mod)
          #t.stat[o,i]<-mod$coefficients[2,3]
          t.stat[o,i]<-mod$coefficients[2,4]
        }
      }
      else {
        #continuous outcome, fit the linear regression
        for(i in 1:ncol(X)){
          data<-data.frame(y=y,x=X[,i])
          mod<-lm(y~x,data=data)
          mod<-summary(mod)
          #t.stat[o,i]<-mod$coefficients[2,3]
          t.stat[o,i]<-mod$coefficients[2,4]
        }
      }
    }
  }
  else {
    for (o in 1:n_O) {
      y=as.vector(Y[,o])
      y.ind=as.vector(Y.ind[,o])
      if (all(is.na(y.ind))) {
        if (length(table(y))==2) {
          #binary outcome, fit the logistic regression
          for(i in 1:ncol(X)){
            data<-data.frame(y=y,x=X[,i])
            mod<-glm(y~x,data=data,family = binomial(link = "logit"))
            mod<-summary(mod)
            #t.stat[o,i]<-mod$coefficients[2,3]
            t.stat[o,i]<-mod$coefficients[2,4]
          }
        }
        else {
          #continuous outcome, fit the linear regression
          for(i in 1:ncol(X)){
            data<-data.frame(y=y,x=X[,i])
            mod<-lm(y~x,data=data)
            mod<-summary(mod)
            #t.stat[o,i]<-mod$coefficients[2,3]
            t.stat[o,i]<-mod$coefficients[2,4]
          }
        }
      } else {
        #survival outcome, fit AFT model with log-logistic parametric family
        for(i in 1:ncol(X)){
          data<-data.frame(y=y,y.ind=y.ind,x=X[,i])
          mod<-survreg(Surv(y, y.ind) ~., data,robust=T,dist="loglogistic")
          mod<-summary(mod)
          #t.stat[o,i]<-mod$coefficients[2,3]
          t.stat[o,i]<-mod$table[2,5]
        }
      }
    }
  }

  if (!is.null(init.NG)) {
    index=apply(t.stat,1,function(x) order(x,decreasing=F)[1:init.NG])
  }
  else {
    if (!is.null(init.pcut)) {
      if (init.p.adjust) {
        index=apply(t.stat,1,function(x) which(p.adjust(x)<init.pcut)) #maybe can do pre-selection on outcomes?
      }
      else {
        index=apply(t.stat,1,function(x) which(x<init.pcut))
      }

    }
    else {
      index=apply(t.stat,1,function(x) order(x,decreasing=F)[1:G/2]) #?? to be carefully think about
    }
  }
  index=unique(unlist(index))
  # print(length(index))
  # if (length(index)==0 | length(index)>(G/2)) {
  #   index=apply(t.stat,1,function(x) order(x,decreasing=F)[1:50])
  #   index=unique(as.vector(index))
  # }
  if (length(index)==0 | length(index)>(G/2)) {
    index=apply(t.stat,1,function(x) order(x,decreasing=F)[1:G/4])
    index=unique(as.vector(index))
  }

  s=rep(0,G)
  s[index]=1 #select the union of top 500 genes highly correlated with either outcome
  #print(sum(s))

  mod.kmeans<-kmeans(X[,index],centers=K,nstart=50) # conduct K-means clustering based on the 500 genes associated with outcome
  C=mod.kmeans$cluster
  Pi=as.vector(prop.table(table(C)))

  #alpha_G=0.9 #given the assumption that the mode-specific clustering are highly consistent with consensus clustering
  alpha_G=alpha_G.init
  alpha_O_vec=alpha_O.init

  L_G=sapply(C,function(c) {
    probs=rep((1-alpha_G)/(K-1),K)
    probs[c]=alpha_G
    return(which(rmultinom(1,1,prob=probs)==1))
  })

  L_O_mat=matrix(NA,nrow=N,ncol=n_O)
  for (o in 1:n_O) {
    alpha_O=alpha_O_vec[o]
    L_O=sapply(C,function(c) {
      probs=rep((1-alpha_O)/(K-1),K)
      probs[c]=alpha_O
      return(which(rmultinom(1,1,prob=probs)==1))
    })
    L_O_mat[,o]=L_O
  }


  MU=sapply(1:K,function(k) {
    center=apply(X[L_G==k,],2,mean)
    #center[s==0]=0 #cluster center based on the L_G clustering and setting the unselected genes (s==0, not associated with outcome) to be zero across clusters
    return(center)
  })

  BETA0=matrix(NA,nrow=n_O,ncol=K)
  BETA=matrix(NA,nrow=n_O,ncol=q-1)

  if (is.null(Y.ind)) {
    for (o in 1:n_O) {
      L_O=as.vector(L_O_mat[,o])
      y=as.vector(Y[,o])
      dat.ini=as.data.frame(cbind(y,L_O,Z[,-1]))
      dat.ini$L_O=factor(dat.ini$L_O,levels=1:K)
      colnames(dat.ini)[1]="y"
      if (length(table(y))==2) {
        #binary outcome, fit the logistic regression
        BETA.temp=as.vector(glm(y~0+.,data=dat.ini,family = binomial(link = "logit"))$coefficients )
        BETA0s=BETA.temp[1:K]
        BETA0[o,]=BETA0s
        BETAs=BETA.temp[-(1:K)]
        BETA[o,]=BETAs
      }
      else {
        #continuous outcome, fit the linear regression
        BETA.temp=as.vector(lm(y~0+.,data=dat.ini)$coefficients )
        BETA0s=BETA.temp[1:K]
        BETA0[o,]=BETA0s
        BETAs=BETA.temp[-(1:K)]
        BETA[o,]=BETAs
      }
    }
  }
  else {
    for (o in 1:n_O) {
      L_O=as.vector(L_O_mat[,o])
      y=as.vector(Y[,o])
      y.ind=as.vector(Y.ind[,o])
      if (all(is.na(y.ind))) {
        dat.ini=as.data.frame(cbind(y,L_O,Z[,-1]))
        dat.ini$L_O=factor(dat.ini$L_O,levels=1:K)
        colnames(dat.ini)[1]="y"
        if (length(table(y))==2) {
          #binary outcome, fit the logistic regression
          BETA.temp=as.vector(glm(y~0+.,data=dat.ini,family = binomial(link = "logit"))$coefficients )
          BETA0s=BETA.temp[1:K]
          BETA0[o,]=BETA0s
          BETAs=BETA.temp[-(1:K)]
          BETA[o,]=BETAs
        }
        else {
          #continuous outcome, fit the linear regression
          BETA.temp=as.vector(lm(y~0+.,data=dat.ini)$coefficients )
          BETA0s=BETA.temp[1:K]
          BETA0[o,]=BETA0s
          BETAs=BETA.temp[-(1:K)]
          BETA[o,]=BETAs
        }
      } else {
        #survival outcome, fit AFT model with log-logistic parametric family
        dat.ini=as.data.frame(cbind(y,y.ind,L_O,Z[,-1]))
        dat.ini$L_O=factor(dat.ini$L_O,levels=1:K)
        colnames(dat.ini)[1:2]=c("y","y.ind")
        BETA.temp<-as.vector(survreg(Surv(y, y.ind) ~0+., dat.ini,robust=T,dist="loglogistic")$coefficients)
        BETA0s=BETA.temp[1:K]
        BETA0[o,]=BETA0s
        BETAs=BETA.temp[-(1:K)]
        BETA[o,]=BETAs
      }
    }
  }

  # w=mean(s)
  # e=c(var(as.vector(MU[s==1,])),var(as.vector(MU[s==0,])))
  tau=apply(X,2,var)

  C.steps=c()
  s.steps=c()
  L_G.steps=c()
  L_O.steps=list()
  BETA.steps=list()
  BETA0.steps=list()
  pars.steps=matrix(NA,nrow=n_iter,ncol=3+n_O+G+n_O+K+1)
  MU.steps=list()
  colnames(pars.steps)=c("alpha_G",paste0("alpha_O_",1:n_O),"e1","e0",paste0("tau",1:G),paste0("sigma_",1:n_O),paste0("pi",1:K),"w")
  for (iter in 1:n_iter) {
    if (iter %% 500==0) print(paste0(iter, "th iteration begin!"))
    updates=Gibbs_step(K=K,MU=MU,BETA0=BETA0,BETA=BETA,s=s,C=C,alpha_G=alpha_G,alpha_O_vec=alpha_O_vec,tau=tau,sigma_vec=sigma_vec,Pi=Pi,L_G,L_O_mat=L_O_mat,
                       beta0_vec=beta0_vec,v0_vec=v0_vec,beta_mat=beta_mat,v_mat=v_mat,X=X,Y=Y,Y.ind=Y.ind,Z=Z,delta=delta,d1=d1,d2=d2,t1=t1,t2=t2,a1=a1,a2=a2,b1=b1,b2=b2,
                       r_G1=r_G1,r_G2=r_G2,r_O1=r_O1,r_O2=r_O2,h=h,cutoff_G = cutoff_G, cutoff_O=cutoff_O)
    MU=updates$MU
    MU.steps=c(MU.steps,list(MU))
    BETA0=updates$BETA0
    BETA0.steps=c(BETA0.steps,list(BETA0))
    BETA.steps=c(BETA.steps,list(updates$BETA))
    s=updates$s
    s.steps=rbind(s.steps,s)
    C=updates$C
    C.steps=rbind(C.steps,C)
    L_G=updates$L_G
    L_G.steps=rbind(L_G.steps,L_G)
    L_O_mat=updates$L_O_mat
    L_O.steps=c(L_O.steps,list(L_O_mat))
    Pi=updates$Pi
    w=updates$w
    sigma_vec=updates$sigma_vec
    tau=updates$tau
    e=updates$e
    alpha_G=updates$alpha_G
    alpha_O_vec=updates$alpha_O_vec
    pars.steps[iter,]=c(alpha_G,alpha_O_vec,e,tau,sigma_vec,Pi,w)
  }

  # saveRDS(list(s.steps=s.steps,C.steps=C.steps,L_G.steps=L_G.steps,L_O.steps=L_O.steps,
  #              BETA0.steps=BETA0.steps,BETA.steps=BETA.steps,MU.steps=MU.steps,pars.steps=pars.steps),
  #         paste0(path,"converge_process.RDS"))
  return(list(s.steps=s.steps,C.steps=C.steps,L_G.steps=L_G.steps,L_O.steps=L_O.steps,
              BETA0.steps=BETA0.steps,BETA.steps=BETA.steps,MU.steps=MU.steps,pars.steps=pars.steps))
  #return(list(res=updates,s.steps=s.steps,C.steps=C.steps,L_G.steps=L_G.steps,L_O.steps=L_O.steps,BETA.steps=BETA.steps))
}

rTbeta=function(n=1,shape1,shape2,cutoff) {
  Ft=pbeta(cutoff,shape1 = shape1,shape2 = shape2)
  Fs=runif(n=n,min=Ft,max=1)
  return(qbeta(Fs,shape1 = shape1,shape2 = shape2))
}
rDir=function(shapes) {
  s_gamma=sapply(shapes,function(x) rgamma(1,shape = x,rate=1))
  return(s_gamma/sum(s_gamma))
}

Gibbs_step=function(K,MU,BETA0,BETA,s,C,alpha_G,alpha_O_vec,tau,sigma_vec,Pi,L_G,L_O_mat,beta0_vec,v0_vec,beta_mat,v_mat,X,Y,Y.ind=NULL,Z,delta,
                       d1,d2,t1,t2,a1,a2,b1,b2,r_G1,r_G2,r_O1,r_O2,h,cutoff_G=1/3,cutoff_O=1/3) {
  G=length(tau)
  N=nrow(Y)
  n_O=ncol(Y)
  q=ncol(Z)

  #given s
  NG=sum(s)
  w=rbeta(1,NG+1,G-NG+1)

  #given MU,s
  e1=1/rgamma(1,shape=NG*K/2+a1,rate=sum(MU[s==1,]^2)/2+a2)
  e0=1/rgamma(1,shape=(G-NG)*K/2+b1,rate=sum(MU[s==0,]^2)/2+b2)
  e=c(e1,e0)
  #given MU,e,w
  sp=(e1)^(-K/2)*exp(rowSums(MU^2)/(-2*e1))
  sq=(e0)^(-K/2)*exp(rowSums(MU^2)/(-2*e0))
  s.update=sapply(1:G,function(g) rbinom(1,1,w*sp[g]/(w*sp[g]+(1-w)*sq[g])))
  # if (sum(s.update)!=0) {s=s.update}
  s=s.update
  #print(sum(s.update))


  #given C,alpha_G,alpha_O,MU,BETA,tau,sigma
  density_LG=t(sapply(1:N,function(i) {
    temp=rep((1-alpha_G)/(K-1),K)
    temp[C[i]]=alpha_G
    return(temp)
  }))

  #density_X=sapply(1:K,function(k) exp(-colSums(((t(X)-MU[,k])^2/tau))/2))*density_LG
  density_X=sapply(1:K,function(k) -colSums(((t(X)-MU[,k])^2/tau))/2)+log(density_LG)
  #density_X=sapply(1:K,function(k) exp(-colSums(((t(X)-MU[,k])^2/tau)[s==1,])/2))*density_LG # only count the density of the selected genes
  prob_LG=density_X-apply(density_X,1,max)
  prob_LG=exp(prob_LG)/rowSums(exp(prob_LG)) #probalematic: some rows are NaN due to density_X is close to 0 for all

  L_G.update=sapply(1:N,function(i) {
    # if (any(is.nan(prob_LG[i,]))) {return(L_G[i])} #not sure whether we solve the issue in this way
    # else {return(which(rmultinom(1,1,prob_LG[i,])==1))}
    # print(i)
    return(which(rmultinom(1,1,prob_LG[i,])==1))
  })
  #if (length(unique(L_G.update))==K) {L_G=L_G.update}
  L_G=L_G.update

  # degenerating=FALSE #debug
  for (o in 1:n_O) {
    alpha_O=alpha_O_vec[o]
    L_O=L_O_mat[,o]
    density_LO=t(sapply(1:N,function(i) {
      temp=rep((1-alpha_O)/(K-1),K)
      temp[C[i]]=alpha_O
      return(temp)
    }))
    BETA0s=as.vector(BETA0[o,])
    BETAs=as.vector(BETA[o,])
    y=as.vector(Y[,o])
    sigma=sigma_vec[o]
    if (!is.null(Y.ind)) {
      y.ind=as.vector(Y.ind[,o])
    }
    else {
      y.ind=rep(NA,length(y))
    }

    if (all(is.na(y.ind))) {
      if (length(table(y))==2) {
        density_y=Z%*%rbind(BETA0s,matrix(rep(BETAs,K),ncol=K))*y-log(1+exp(Z%*%rbind(BETA0s,matrix(rep(BETAs,K),ncol=K))))
        density_y=exp(density_y)*density_LO
      }
      else {
        density_y=exp((Z%*%rbind(BETA0s,matrix(rep(BETAs,K),ncol=K))-y)^2/(-2*sigma))*density_LO
      }
    }
    else {
      # density_y=exp(1/sqrt(sigma)*(y-Z%*%rbind(BETA0s,matrix(rep(BETAs,K),ncol=K))))/(1+exp(1/sqrt(sigma)*(y-Z%*%rbind(BETA0s,matrix(rep(BETAs,K),ncol=K)))))^(1+y.ind)*density_LO
      density_y=1/sqrt(sigma)*(y-Z%*%rbind(BETA0s,matrix(rep(BETAs,K),ncol=K)))-log(1+exp(1/sqrt(sigma)*(y-Z%*%rbind(BETA0s,matrix(rep(BETAs,K),ncol=K)))))*(1+y.ind)
      density_y=exp(density_y)*density_LO
    }

    prob_LO=density_y/rowSums(density_y)

    #debug
    # if (any(is.nan(prob_LO))) {
    #   degenerating=TRUE
    #   break
    # }
    #

    L_O.update=sapply(1:N,function(i) {
      return(which(rmultinom(1,1,prob_LO[i,])==1))
    })
    L_O=L_O.update
    L_O_mat[,o]=L_O
    alpha_O_vec[o]=rTbeta(n=1,shape1=r_O1[o]+sum(L_O==C),shape2=r_O2[o]+sum(L_O!=C),cutoff=cutoff_O[o])
  }

  #debug
  # if (degenerating) {
  #   return(list(density_y=density_y,density_LO=density_LO,
  #               fy=exp((Z%*%rbind(BETA0s,matrix(rep(BETAs,K),ncol=K))-y)^2/(-2*sigma)),
  #               BETA0s=BETA0s,
  #               BETAs=BETAs,
  #               sigma=sigma))
  # }
  #
  #given L_O,L_G
  alpha_G=rTbeta(n=1,shape1=r_G1+sum(L_G==C),shape2=r_G2+sum(L_G!=C),cutoff=cutoff_G)

  #given L_O,L_G,alpha_G,alpha_O,Pi
  prob_LG=sapply(1:N,function(i) {
    temp=rep((1-alpha_G)/(K-1),K)
    temp[L_G[i]]=alpha_G
    return(temp)
  })

  alpha_O=alpha_O_vec[1]
  L_O=L_O_mat[,1]
  prob_LO=sapply(1:N,function(i) {
    temp=rep((1-alpha_O)/(K-1),K)
    temp[L_O[i]]=alpha_O
    return(temp)
  })
  prob_C=prob_LG*prob_LO

  if (n_O>1) {
    for (o in 2:n_O) {
      alpha_O=alpha_O_vec[o]
      L_O=L_O_mat[,o]
      prob_LO.temp=sapply(1:N,function(i) {
        temp=rep((1-alpha_O)/(K-1),K)
        temp[L_O[i]]=alpha_O
        return(temp)
      })
      prob_C=prob_C*prob_LO.temp
    }
  }

  prob_C=t(prob_C*Pi)
  C.update=sapply(1:N,function(i) {
    which(rmultinom(1,1,prob_C[i,]/sum(prob_C[i,]))==1)
  })

  #if (length(unique(C.update))==K) {C=C.update}
  C=C.update

  #given C
  C_n=sapply(1:K,function(k) sum(C==k))
  Pi=rDir(shapes = h+C_n)

  #given s,e,L_G,tau
  es=e[2-s]
  MU=do.call(cbind,lapply(1:K,function(k) {
    sd_Vec=1/(1/es+sum(L_G==k)*(1/tau))
    m_Vec=diag(sd_Vec*1/tau) %*% colSums(X[L_G==k,,drop=F])
    mu_k=mvrnorm(mu=m_Vec,Sigma=diag(sd_Vec))
    return(mu_k)
  }))

  #given L_O,sigma
  # infBeta0=FALSE #debug
  if (is.null(Y.ind)) {
    for (o in 1:n_O) {
      alpha_O=alpha_O_vec[o]
      L_O=L_O_mat[,o]
      v0=v0_vec[o] #prior parameters,will not be updated?
      beta0=beta0_vec[o] #prior parameters,will not be updated?
      sigma=sigma_vec[o] #prior parameters,will not be updated?
      y=as.vector(Y[,o])
      BETAs=as.vector(BETA[o,])
      if (length(table(y))==2) {
        #binary outcome
        omega=sapply(1:length(L_O),function(i) BayesLogit::rpg(num=1, h=1, z=BETA0[o,L_O[i]]+Z[i,-1,drop=F]%*%BETAs)) #generate Polya Gamma variables
        BETA0s=sapply(1:K,function(k) {
          if (sum(L_O==k)==0) {
            beta_k=beta0
          }
          else {
            sd0=1/(1/(sigma*v0)+sum(omega[L_O==k]))
            m0=sd0*(1/(sigma*v0)*beta0+sum(y[L_O==k]-Z[L_O==k,-1,drop=F]%*%BETAs*omega[L_O==k]-1/2))
            beta_k=rnorm(1,mean=m0,sd=sd0)
          }
          return(beta_k)
        })
        BETA0[o,]=BETA0s

        v_vec=as.vector(v_mat[o,])
        beta_vec=as.vector(beta_mat[o,])

        sd_Mat=solve(diag(1/v_vec)*(1/sigma)+t(Z[,-1,drop=F])%*%diag(omega)%*%Z[,-1,drop=F])
        m_Vec=sd_Mat%*% (1/v_vec*(1/sigma)*beta_vec+t(Z[,-1,drop=F])%*%(y-1/2-BETA0s[L_O]*omega))
        BETAs=mvrnorm(mu=m_Vec,Sigma=sd_Mat)
        BETA[o,]=BETAs

        #given L_O,BETA
        sigma=1/rgamma(n=1,shape=t1+K/2+(q-1)/2,rate=t2+sum((BETAs-beta_vec)^2/v_vec)/2+sum((BETA0s-beta0)^2/v0)/2)
        sigma_vec[o]=sigma
      }
      else {
        #continuous outcome
        BETA0s=sapply(1:K,function(k) {
          if (sum(L_O==k)==0) {
            beta_k=beta0
          }
          else {
            sd0=v0/(1+v0*sum(L_O==k))
            m0=sd0*sum(y[L_O==k]-Z[L_O==k,-1,drop=F]%*%BETAs)+beta0/(1+v0*sum(L_O==k))
            beta_k=rnorm(1,mean=m0,sd=sd0*sigma)
          }
          return(beta_k)
        })
        BETA0[o,]=BETA0s

        v_vec=as.vector(v_mat[o,])
        beta_vec=as.vector(beta_mat[o,])

        sd_Mat=solve(diag(1/v_vec)+t(Z[,-1,drop=F])%*%Z[,-1,drop=F])
        m_Vec=sd_Mat%*% (1/v_vec*beta_vec+t(Z[,-1,drop=F])%*%(y-BETA0s[L_O]))
        BETAs=mvrnorm(mu=m_Vec,Sigma=sd_Mat*sigma)
        BETA[o,]=BETAs

        #given L_O,BETA
        sigma=1/rgamma(n=1,shape=t1+N/2+K/2+(q-1)/2,rate=t2+sum((y-BETA0s[L_O]-Z[,-1,drop=F]%*%BETAs)^2)/2+sum((BETAs-beta_vec)^2/v_vec)/2+sum((BETA0s-beta0)^2/v0)/2)
        sigma_vec[o]=sigma
      }
    }
  }
  else {
    for (o in 1:n_O) {
      alpha_O=alpha_O_vec[o]
      L_O=L_O_mat[,o]
      v0=v0_vec[o]
      beta0=beta0_vec[o]
      sigma=sigma_vec[o]
      y=as.vector(Y[,o])
      y.ind=as.vector(Y.ind[,o])
      BETAs=as.vector(BETA[o,])
      if (all(is.na(y.ind))) {
        if (length(table(y))==2) {
          #binary outcome
          omega=sapply(1:length(L_O),function(i) BayesLogit::rpg(num=1, h=1, z=BETA0[o,L_O[i]]+Z[i,-1,drop=F]%*%BETAs)) #generate Polya Gamma variables
          BETA0s=sapply(1:K,function(k) {
            if (sum(L_O==k)==0) {
              beta_k=beta0
              # beta_k=BETA0[o,k]
            }
            else {
              sd0=1/(1/(sigma*v0)+sum(omega[L_O==k]))
              m0=sd0*(1/(sigma*v0)*beta0+sum(y[L_O==k]-Z[L_O==k,-1,drop=F]%*%BETAs*omega[L_O==k]-1/2))
              beta_k=rnorm(1,mean=m0,sd=sd0)
            }
            return(beta_k)
          })
          BETA0[o,]=BETA0s

          v_vec=as.vector(v_mat[o,])
          beta_vec=as.vector(beta_mat[o,])

          sd_Mat=solve(diag(1/v_vec)*(1/sigma)+t(Z[,-1,drop=F])%*%diag(omega)%*%Z[,-1,drop=F])
          m_Vec=sd_Mat%*% (1/(v_vec*sigma)*beta_vec+t(Z[,-1,drop=F])%*%(y-1/2-BETA0s[L_O]*omega))
          BETAs=mvrnorm(mu=m_Vec,Sigma=sd_Mat)
          BETA[o,]=BETAs

          #given L_O,BETA
          sigma=1/rgamma(n=1,shape=t1+K/2+(q-1)/2,rate=t2+sum((BETAs-beta_vec)^2/v_vec)/2+sum((BETA0s-beta0)^2/v0)/2)
          sigma_vec[o]=sigma
        }
        else {
          #continuous outcome
          BETA0s=sapply(1:K,function(k) {
            if (sum(L_O==k)==0) {
              beta_k=beta0
              # beta_k=BETA0[o,k]
            }
            else {
              sd0=v0/(1+v0*sum(L_O==k))
              m0=sd0*sum(y[L_O==k]-Z[L_O==k,-1,drop=F]%*%BETAs)+beta0/(1+v0*sum(L_O==k))
              beta_k=rnorm(1,mean=m0,sd=sd0*sigma)
            }
            return(beta_k)
          })
          BETA0[o,]=BETA0s

          v_vec=as.vector(v_mat[o,])
          beta_vec=as.vector(beta_mat[o,])

          sd_Mat=solve(diag(1/v_vec)+t(Z[,-1,drop=F])%*%Z[,-1,drop=F])
          m_Vec=sd_Mat%*% (1/v_vec*beta_vec+t(Z[,-1,drop=F])%*%(y-BETA0s[L_O]))
          BETAs=mvrnorm(mu=m_Vec,Sigma=sd_Mat*sigma)
          BETA[o,]=BETAs

          #given L_O,BETA
          sigma=1/rgamma(n=1,shape=t1+N/2+K/2+(q-1)/2,rate=t2+sum((y-BETA0s[L_O]-Z[,-1,drop=F]%*%BETAs)^2)/2+sum((BETAs-beta_vec)^2/v_vec)/2+sum((BETA0s-beta0)^2/v0)/2)
          sigma_vec[o]=sigma
        }
      }
      else {
        #survival outcome
        omega=sapply(1:length(L_O),function(i) BayesLogit::rpg(num=1, h=1+y.ind[i], z=1/sqrt(sigma)*(log(y[i])-BETA0[o,L_O[i]]-Z[i,-1,drop=F]%*%BETAs))) #generate Polya Gamma variables
        BETA0s=BETA0[o,]
        for (k in 1:K) {
          if (sum(L_O==k)==0) {
            BETA0s[k]=rnorm(1,mean=beta0,sd=v0*sigma)
          }
          else {
            sd0=sigma/(1/v0+sum(omega[L_O==k]))
            m0=sd0*(beta0/(v0*sigma)+sum((log(y)[L_O==k]-Z[L_O==k,-1,drop=F]%*%BETAs)*omega[L_O==k]/sigma-(1-y.ind[L_O==k])/(2*sqrt(sigma))))
            BETA0s[k]=rnorm(1,mean=m0,sd=sd0)
          }
        }
        BETA0[o,]=BETA0s

        v_vec=as.vector(v_mat[o,])
        beta_vec=as.vector(beta_mat[o,])

        sd_Mat=sigma*solve(diag(1/v_vec)+t(Z[,-1,drop=F])%*%diag(omega)%*%Z[,-1,drop=F])
        m_Vec=sd_Mat%*% (1/(v_vec*sigma)*beta_vec+t(Z[,-1,drop=F])%*%((log(y)-BETA0s[L_O])*omega/sigma-(1-y.ind)/(2*sqrt(sigma))))
        BETAs=mvrnorm(mu=m_Vec,Sigma=sd_Mat)
        BETA[o,]=BETAs

        #given L_O,BETA to update sigma(Metropolis-Hasting sampler)
        sigma_condidate=(exp(log(sqrt(sigma))+delta*rnorm(1,mean=0,sd=1)))^2
        BETA0lr=sum(log(dnorm(BETA0s,mean=beta0,sd=v0*sigma_condidate)))-
          sum(log(dnorm(BETA0s,mean=beta0,sd=v0*sigma)))
        omegalr=sum(log(1+exp((log(y)-BETA0s[L_O]-Z[,-1,drop=F]%*%BETAs)/sqrt(sigma_condidate)))-
                      log(1+exp((log(y)-BETA0s[L_O]-Z[,-1,drop=F]%*%BETAs)/sqrt(sigma)))+
                      (1/sqrt(sigma)-1/sqrt(sigma_condidate))/2*(log(y)-BETA0s[L_O]-Z[,-1,drop=F]%*%BETAs)+
                      (1/sigma-1/sigma_condidate)/2*(log(y)-BETA0s[L_O]-Z[,-1,drop=F]%*%BETAs)^2*omega)
        BETAlr=log(dmvnorm(BETA[o,,drop=F],beta_vec,diag(v_vec)*sigma_condidate))-log(dmvnorm(BETA[o,,drop=F],beta_vec,diag(v_vec)*sigma))
        ylr=sum(y.ind)*(log(sqrt(sigma))-log(sqrt(sigma_condidate)))+
          sum((1/sqrt(sigma_condidate)-1/sqrt(sigma))*(log(y)-BETA0s[L_O]-Z[,-1,drop=F]%*%BETAs)+
                (1+y.ind)*(log(1+exp((log(y)-BETA0s[L_O]-Z[,-1,drop=F]%*%BETAs)/sqrt(sigma)))-
                             log(1+exp((log(y)-BETA0s[L_O]-Z[,-1,drop=F]%*%BETAs)/sqrt(sigma_condidate)))))
        llr=BETA0lr+BETAlr+omegalr+ylr+
          log(dgamma(1/sigma_condidate,shape=t1,rate=t2))-log(dgamma(1/sigma,shape=t1,rate=t2))
        p_update=exp(llr)*(sqrt(sigma_condidate)/sqrt(sigma))
        if (runif(1)<=p_update) {
          sigma_vec[o]=sigma_condidate
        }

      }

    }
  }

  #debug
  # if (infBeta0) {
  #   return(list(o=o,L_O_mat=L_O_mat,v0_vec=v0_vec,
  #               sigma_vec=sigma_vec,beta0_vec=beta0_vec,BETA=BETA))
  # }
  #

  #given L_G,MU
  tau=sapply(1:G,function(g) {
    1/(rgamma(n=1,shape=d1+N/2,rate=d2+sum((X[,g]-MU[g,L_G])^2)/2))
  })

  return(list(MU=MU,BETA0=BETA0,BETA=BETA,s=s,C=C,L_G=L_G,L_O_mat=L_O_mat,Pi=Pi,w=w,sigma_vec=sigma_vec,tau=tau,e=e,alpha_G=alpha_G,alpha_O_vec=alpha_O_vec))
}

