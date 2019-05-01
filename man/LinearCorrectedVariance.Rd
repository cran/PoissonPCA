\name{LinearCorrectedVariance}
\alias{LinearCorrectedVariance}
\alias{LinearCorrectedVarianceSeqDepth}
\alias{TransformedVariance}
\alias{TransformedVarianceECV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimates variance of a transformation of latent Poisson means
%%  ~~function to do ... ~~
}
\description{
  Given a data matrix X[i,j] which follows a Poisson distribution with
  mean Lambda[i,j], this function estimates the covariance matrix of a
  transformation f of the latent Lambda.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LinearCorrectedVariance(X)
LinearCorrectedVarianceSeqDepth(X)
TransformedVariance(X,g,CVar)
TransformedVarianceECV(X,g,ECVar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{the data matrix}

  \item{g}{an estimator of the transformation function f(Lambda). That
    is, if X~Poisson(Lambda), g(X) should be an estimator of f(Lambda).}

  \item{CVar}{an estimator of the conditional variance of g(X)
    conditional on Lambda.}

  \item{ECVar}{an estimator of the conditional variance of g(X)
    conditional on Lambda.}
}
\details{
\code{LinearCorrectedVariance} merely estimates the covariance matrix of
    the latent Poisson means without
    transformation. \code{LinearCorrectedVarianceSeqDepth} deals with
    the common case in microbiome and other analysis, where the Poisson
    means are subject to large multiplicative noise not associated with
    the parameters of interest. In these cases, we would like to
    estimate the covariance of the compositional form of Lambda. That
    is, we want to scale the rows of Lambda to all have sum 1, and
    estimate the covariance matrix of the resultant matrix. This method
    uses the actual row sums of X as estimates of the scaling to be
    performed. \code{TransformedVariance} estimates the variance of
    a function of Lambda. It takes two additional parameters: \code{g}
    and \code{CVar} which are functions of X. \code{g} should be an
    estimator for the desired transformation f(Lambda) from an
    observation X. For example, if f(Lambda)=Lambda^2, then the
    unbiassed estimator is X*(X-1). \code{CVar} is an estimator for the
    conditional variance of g(X) given Lambda. For example, if
    f(Lambda)=Lambda^2, and we use the unbiassed g(X)=X*(X-1), then the
    variance of g(X) is 4*Lambda^3+3*Lambda^2, so an unbiassed estimator
    for this is
    \code{CVar(X)}=4*X*(X-1)*(X-2)+3*X*(X-1)=X*(X-1)*(4*X-5). The
    function \code{polynomial_transformation} will compute the unbiassed
    estimators for a given polynomial. The
    function \code{makelogtransformation} compute estimators for the log
    function. \code{TransformedVarianceECV} is the same as
    \code{TransformedVariance}, except that the third parameter
    estimates the average conditional variance from a sample of values
    of X, rather than a single value. 
}
\value{
  An estimated covariance matrix for the transformed latent means.
}
%\references{}
\author{Toby Kenney \email{tkenney@mathstat.dal.ca} and Tianshu Huang \email{}}
\examples{

n<-20  #20 observations
p<-5   #5 dimensions
r<-2   #rank 2

mean<-10*c(1,3,2,1,1)

Z<-rnorm(n*r)
dim(Z)<-c(n,r)
U<-rnorm(p*r)
dim(U)<-c(r,p)

Latent<-Z\%*\%U+rep(1,n)\%*\%t(mean)

X<-rpois(n*p,as.vector(Latent))
dim(X)<-c(n,p)

LinearCorrectedVariance(X)

seqdepth<-exp(rnorm(n)+2)
Xseqdep<-rpois(n*p,as.vector(diag(seqdepth)\%*\%Latent))
dim(Xseqdep)<-c(n,p)

LinearCorrectedVarianceSeqDepth(Xseqdep)


squaretransform<-polynomial_transformation(c(1,0))

Xsq<-rpois(n*p,as.vector(diag(seqdepth)\%*\%Latent)^2)

dim(Xsq)<-c(n,p)

TransformedVariance(Xsq,squaretransform$g,squaretransform$CVar)

Xexp<-rpois(n*p,as.vector(diag(seqdepth)\%*\%exp(Latent)))

logtrans<-makelogtransformation(3,4)

TransformedVarianceECV(X,logtrans$g,logtrans$ECVar)


}
\keyword{PCA}
