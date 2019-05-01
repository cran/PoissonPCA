\name{Poisson_Corrected_PCA}
\alias{Poisson_Corrected_PCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PCA with Poisson measurement error
  %% ~~function to do ... ~~
}
\description{Estimates the principal
components of the latent Poisson means (possibly with transformation) of
high-dimensional data with independent Poisson measurement error.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Poisson_Corrected_PCA(X,k=dim(X)[2]-1,transformation=NULL,seqdepth=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Matrix or data frame of count variables}
  \item{k}{Number of principal components to calculate.}
  \item{transformation}{For estimating the principal components of a
    transformation of the Poisson mean.}
  \item{seqdepth}{Indicates what sort of sequencing depth correction
    should be used (if any).}
}
\details{
   The options for the transformation parameter are:

    NULL or "linear" - these perform no transformation.

    "log" - this performs a logarithmic transformation

    a list of the following functions:

        f(x) - evaluates the function
	deriv(x) - evaluates the derivative of the function
	solvefunction(target) - evaluates the inverse of the function
	g(x) - an estimator for f(lambda) from a Poisson observation x
	          with mean lambda
	CVar(x) - an estimator for the conditional variance of g(x)
	          conditional on lambda from the observed value x

    the function polynomial_transformation creates such a list in the
    case where f is a polynomial using unbiassed estimators for g and
    CVar. The function makelogtransformation creates an estimator for
    the logarithmic transformation. The "log" option uses this function
    with parameters a=3 and N=4, which from experiments appear to
    produce reasonable results in most situations.
		  

    The options for the seqdepth parameter are:

    FALSE - indicating no sequencing depth correction

    TRUE - indicating standard sequencing depth correction for linear
             PCA

    "minvar" - uses the minimum covariance estimator for the corrected
              variance. This subtracts the largest constant from all
              entries of the matrix, such that the matrix is still
              non-negative definite.	     

    "compositional" - uses the best compositional variance approximation
                      to the estimated covariance matrix.


    The package estimates latent principal components using the methods
    in http://arxiv.org/abs/1904.11745
}


\value{
  An object of type "princomp" and "transformedprincomp" that has the
  following components:
  \item{sdev}{The standard deviation associated to each principal component}
  \item{loadings}{The principal component vectors}
  \item{center}{The mean of the transformed data}
  \item{scale}{A vector of ones of length n}
  \item{n.obs}{The number of observations}
  \item{scores}{The principal scores. For the linear transformation,
  these are just the projection of the data onto the principal component
  space. For transformed principal components, these use a combination
  of likelihood and mean squared error.}
  \item{means}{The corresponding estimated untransformed Poisson
  means. This provides a useful diagnostic of the performance in
  simulation studies. These means should be closer to the true Lambda
  than the original X data.}
  \item{variance}{The corrected covariance matrix for the transformed
    latent Sigma.}
  \item{non_compositional_variance}{The corrected covariance matrix
    without sequencing depth correction.}
  \item{call}{The function call used}
}
%\references{}
\author{Toby Kenney \email{tkenney@mathstat.dal.ca} and Tianshu Huang \email{}}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{}
\examples{

set.seed(12345)
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

Poisson_Corrected_PCA(X,k=2,transformation=NULL,seqdepth=FALSE)

seqdepth<-exp(rnorm(n)+2)
Xseqdep<-rpois(n*p,as.vector(diag(seqdepth)\%*\%Latent))
dim(Xseqdep)<-c(n,p)

Poisson_Corrected_PCA(Xseqdep,k=2,transformation=NULL,seqdepth=TRUE)

squaretransform<-polynomial_transformation(c(1,0))

Xexp<-rpois(n*p,as.vector(diag(seqdepth)\%*\%exp(Latent)))

Poisson_Corrected_PCA(Xseqdep,k=2,transformation="log",seqdepth="minvar")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{PCA}
\keyword{Poisson measurement error}% __ONLY ONE__ keyword per line

