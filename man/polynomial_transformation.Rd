\name{polynomial_transformation}
\alias{polynomial_transformation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{constructs a polynomial transformation for use with functions
  from the PoissoncorrectedPCA package.
%%  ~~function to do ... ~~
}
\description{When we are dealing with a transformation of the latent
  Poisson mean Lambda, we need various useful functions. This function
  computes the necessary functions for a polynomial, and returns a list
  of the required functions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
polynomial_transformation(coeffs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coeffs}{A vector of coefficents of the polynomial. The constant
    term should not be included.}
}
\details{
   The coefficients of the polynomial should be given in order of
    decreasing degree, and should not include the constant term. For
    example \code{"coeffs"=c(1,2,3)} refers to the polynomial
    X^3+2*X^2+3*X. This function returns a list of functions for dealing
    with this transformation.
}
\value{
  \item{f}{evaluates the transformation}
  \item{g}{an estimator for the transformation of a latent Poisson mean}
  \item{solve}{computes the inverse transformation (often used for simulations)}
  \item{CVar}{an estimator for the conditional variance of g(X)}
}
%\references{}
\author{Toby Kenney \email{tkenney@mathstat.dal.ca} and Tianshu Huang \email{}}
\examples{

cubic<-polynomial_transformation(c(1,0,0))

X<-rpois(100,1.8^3)

gX<-X
varX<-X
for(i in 1:100){
gX[i]<-cubic$g(X[i])
varX[i]<-cubic$CVar(X)
}
mean(gX)
var(gX)
mean(varX)



}
\keyword{PCA}


