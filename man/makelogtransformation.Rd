\name{makelogtransformation}
\alias{makelogtransformation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{constructs a log transformation for use with functions
  from the PoissoncorrectedPCA package.
%%  ~~function to do ... ~~
}
\description{When we are dealing with a transformation of the latent
  Poisson mean Lambda, we need various useful functions. This function
  computes the necessary functions for the log transformation, and returns a list of the required functions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makelogtransformation(a,N,uselog=6,unbiassed=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{The point about which to expand the Taylor series (see details)}
  \item{N}{The number of terms in the Taylor series to expand (see
    details)}
  \item{uselog}{Value above which we use the logarithm to approximate
    g(x). Typically should not be larger than 2a.}
  \item{unbiassed}{Indicates that the recommended unbiassed method
    should be used.}
}
\details{
  The logarithmic transformation is fundamentally unestimable. There is
  no estimator which is an unbiassed estimator for log(Lambda). This is
  because the logarithm function has a singularity at zero, so has no
  globally convergent Taylor series expansion. Instead, we aim to use an
  approximately unbiassed estimator. For large enough X, g(X)=log(X) is
  a reasonable estimator. For smaller X, we need to compute a Taylor
  series for exp(-Lambda)log(Lambda). We do this from the equation
  log(x)=log(a)+log(x/a) and the Taylor expansion
  log(1+y)=y-y^2/2+y^3/3-... where y=x/a-1. This has radius of
  convergence 1, so will converge provided 0<x<2a. However, if we try to
  convert it to a polynomial in x, the coefficients will
  diverge. Instead, we truncate this Taylor series in y at a chosen
  number N terms. If the x is close to a, this truncated Taylor series
  should give an approximately unbiassed estimator for log(Lambda).
  Choice of N can have some effect. Larger values of N reduce the bias
  of g(X) but increase the variance. Experimentally, a=3 and N=6 seem to
  produce reasonable results, with g(X)=log(X) for X>6.
}
\value{
  \item{type}{="log"}
  \item{f}{function which evaluates the transformation}
  \item{g}{an estimator for the transformation of a latent Poisson mean}
  \item{solve}{function which computes the inverse transformation (often used for simulations)}
  \item{ECVar}{an estimator for the average conditional variance of g(X)}
}
%\references{}
\author{Toby Kenney \email{tkenney@mathstat.dal.ca} and Tianshu Huang}
\examples{

logtrans<-makelogtransformation(5,6)
X<-rpois(100,exp(1.4))
gX<-X
for(i in 1:100){
gX[i]<-logtrans$g(X[i])
}
mean(gX)
var(gX)
logtrans$ECVar(X)

}
\keyword{PCA}


