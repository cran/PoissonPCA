\name{get_scores}
\alias{get_scores}
\alias{get_scores_log}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates principal scores for Poisson-noise corrected PCA
%%  ~~function to do ... ~~
}
\description{
  This function is based on principal component analysis of a
  transformation of latent Poisson means of a sample. Given the
  estimated principal components of the latent Poisson means, this
  function estimates scores using a combination of likelihood and mean
  squared error.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_scores(X,V,d,k,transformation,mu)
get_scores_log(X,V,d,k,mu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The data matrix}
  \item{V}{Vector of all principal components of the transformed latent means}
  \item{d}{Eigenvalues corresponding to the principal components}
  \item{k}{Number of principal components to project onto}
  \item{transformation}{The transformation to be applied to the latent means}
  \item{mu}{The mean of the transformed latent means}
}
\details{
  This function estimates the latent transformed Poisson means in order
  to minimise a combination of the log-likelihood plus the squared
  residuals of the projection of these latent means onto the first k
  principal components. Note that for transformed Poisson PCA, the
  scores are not nested, so the choice of k will have an impact on the
  projection. The \code{get_scores_log} function deals with the special
  case where the transformation is the log function.
}
\value{
  \item{scores}{The principal scores}
  \item{means}{The corresponding estimated latent Poisson means}
}
%\references{}
\author{Toby Kenney \email{tkenney@mathstat.dal.ca} and Tianshu Huang \email{}}
\examples{

n<-20  #20 observations
p<-5   #5 dimensions
r<-2   #rank 2

mean<-10*c(1,3,2,1,1)

set.seed(12345)

Z<-rnorm(n*r)
dim(Z)<-c(n,r)
U<-rnorm(p*r)
dim(U)<-c(r,p)

Latent<-Z\%*\%U+rep(1,n)\%*\%t(mean)

X<-rpois(n*p,as.vector(Latent))
dim(X)<-c(n,p)

Sigma<-LinearCorrectedVariance(X[-n,])

eig<-eigen(Sigma)


get_scores(X[n,],eig$vectors,eig$values,r,"linear",colMeans(X[-n,]))

Xlog<-rpois(n*p,exp(as.vector(Latent)+3))
dim(Xlog)<-c(n,p)

logtrans<-makelogtransformation(3,4)

Sigmalog<-TransformedVarianceECV(X[-n,],logtrans$g,logtrans$ECVar)

eiglog<-eigen(Sigmalog)

gX<-X[-n,]

if(!is.null(logtrans)){
    for(i in 1:(n-1)){
        for(j in 1:p){
            gX[i,j]<-logtrans$g(X[i,j])
        }
    }
}                    

mu<-colMeans(gX)        


get_scores_log(X[n,],eiglog$vectors,eiglog$values,r,mu)

}
\keyword{PCA}



