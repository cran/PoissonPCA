\name{make_compositional_variance}
\alias{make_compositional_variance}
\alias{make_compositional_min_var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Converts a covariance matrix to compositional form
%%  ~~function to do ... ~~
}
\description{Given a covariance matrix, removes multiplicative noise 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_compositional_variance(Sigma)
make_compositional_min_var(Sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigma}{the uncorrected covariance matrix}
}
\details{
  The two functions use different
  methods. \code{make_compositional_variance} calculates the variance of
  compositional data that agrees with Sigma (viewed as a bilinear form)
  on compositional vectors. That is, the return value Sigma_c is a
  symmetric matrix which satisfies \code{t(u)\%*\%Sigma_c\%*\%v=t(u)\%*\%Sigma\%*\%v} 
  for any compositional vectors u and v, and also \code{rowSums(Sigma_c)=0}.
}
\value{
  The compositionally corrected covariance matrix.
}
%\references{}
\author{Toby Kenney \email{tkenney@mathstat.dal.ca} and Tianshu Huang \email{}}
\examples{

n<-10
p<-5

X<-rnorm(n*p)

dim(X)<-c(n,p)
Sigma<-t(X)\%*\%X/(n-1)

SigmaComp<-make_compositional_variance(Sigma)

SigmaCompMin<-make_compositional_min_var(Sigma)

}
\keyword{PCA}


