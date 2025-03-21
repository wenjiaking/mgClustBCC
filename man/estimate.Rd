% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.R
\name{estimate}
\alias{estimate}
\title{Point estimates based on posterior samples}
\usage{
estimate(MCMCres, nsamples = 1000, K)
}
\arguments{
\item{MCMCres}{output of MCMC function}

\item{nsamples}{the number of posterior samples after burn-in (must smaller than n_iter in MCMC function), 1000 by default}

\item{K}{the number of clusters}
}
\value{
A list of point estimates.
\itemize{
\item{s.est: }{a vector of length G, indicating whether the feature is selected (1) or not (0)}
\item{alpha.est: }{a vector contains the estiamtes of omics-specific consensus rate and outcome-specific consensus rates}
\item{pi.est: }{a vector of length K, the estiamtes of cluster probabilities}
\item{MU.est: }{a G by K matrix where each row represents the omics-specific cluster centers for each feature}
\item{BETA0.est: }{an O by K matrix where each row represents the intercepts across outcome-specific cluster for each outcome}
\item{BETA.est: }{an O by q matrix where each row represents the covariate coefficients for each outcome}
\item{C.est: }{a vector of length N, the master cluster estimate for each sample}
\item{L_G.est: }{a vector of length N, the omics-specific cluster estimate for each sample}
\item{L_O.est: }{a N by O matrix where each row represents the outcome-specific cluster estimates for each sample}
\item{TAU.est: }{a vector of length G, the estimates of omics variance}
\item{sigma.est: }{a vector of length O, the estimates of outcome variance}
\item{e.est: }{a vector of length 2, the variance estiamtes of cluster-varying and stable feature}
}
}
\description{
Point estimates based on posterior samples
}
\examples{
\dontrun{
  test.est=estimate(test,nsamples = 500,K=3)
}
}
