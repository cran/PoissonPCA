
squarepolynomial<-function(coeffs){
    degree<-length(coeffs)
    squarecoeffs<-rep(0,2*degree)
    for(i in seq_len(degree)){
        squarecoeffs[i]<-sum(coeffs[1:i]*coeffs[i:1])
    }
    for(i in seq_len(degree)){
        squarecoeffs[i+degree-1]<-sum(coeffs[i:degree]*coeffs[degree:i])
    }
    return(squarecoeffs)
}

polynomial_transformation<-function(coeffs){
#coeff is a vector of coefficients from highest to lowest.
#constant coefficient is assumed to be 0    
    force(coeffs)
    degree<-length(coeffs)
    if(degree>1){
        increasingcoeffs<-coeffs[degree:1]
        f<-function(x){#x needs to be a scalar
            return(sum((x^(degree:1))*coeffs))
        }
        deriv<-function(x){#x needs to be a scalar
            return(sum((x^((degree-1):1))*coeffs[1:(degree-1)]*(degree:2))+coeffs[degree])
        }
        solvefunction<-function(target){
            x<-(target/coeffs[1])^(1/degree)
            val<-sum((x^(degree:1))*coeffs)-target
            while(val*val>1e-16){
                der<-(sum((x^((degree-1):1))*coeffs[1:(degree-1)]*(degree:2))+coeffs[degree])
                x<-x-val/der
                val<-sum((x^(degree:1))*coeffs)-target            
            }
            return(x)
        }    
        g<-function(x){
            answer<-x
            for(i in 1:length(x)){
                powerests<-cumprod(x[i]-(0:(degree-1)))
                answer[i]<-sum(increasingcoeffs*powerests)
            }
            return(answer)
        }
        squarecoeffs<-squarepolynomial(coeffs)
        squarecoeffs<-c(0,squarecoeffs[(2*degree-1):1])
        CVar<-function(x){
            powerests<-cumprod(x-(0:(2*degree-1)))
            return((sum(increasingcoeffs*powerests[1:degree]))^2-sum(squarecoeffs*powerests))
        }
    }else if(degree==1){
        ## Linear polynomial

        f<-function(x){#x needs to be a scalar
            return(x*coeffs)
        }
        deriv<-function(x){#x needs to be a scalar
            return(coeffs)
        }
        solvefunction<-function(target){
            return(target/coeffs)
        }    
        g<-function(x){
            return(x)
        }
        CVar<-function(x){
            return(coeffs^2*x*(x-1))
        }        
    }else{
        return(NULL)
    }
    return(list("f"=f,"g"=g,"solve"=solvefunction,"CVar"=CVar,"type"="polynomial"))
}


TransformedVariance<-function(X,g,CVar){
    #Estimates the variance of f(lambda)
    #X is the data matrix
    #g(x) is an estimator of the transformed Poisson mean
    #CVar(lambda) is the expected conditional variance of g(X) given lambda
    #This does not include any sequencing-depth correction
    
    gX<-rep(0,(dim(X)[1])*(dim(X)[2]))
    dim(gX)<-dim(X)
    CVX<-rep(0,(dim(X)[1])*(dim(X)[2]))
    dim(CVX)<-dim(X)
    for(i in seq_len(dim(X)[1])){
        for(j in seq_len(dim(X)[2])){
            gX[i,j]=g(X[i,j])
            CVX[i,j]=CVar(X[i,j])
        }
    }
    vargx<-stats::var(gX)
    
    ECVar<-colMeans(CVX)
    varf<-stats::var(gX)-diag(ECVar)
    return(varf)
}



TransformedVarianceECV<-function(X,g,ECVar){
    #Estimates the variance of f(lambda)
    #X is the data matrix
    #g(x) is an estimator of the transformed Poisson mean
    #CVar(lambda) is the expected conditional variance of g(X) given lambda
    #This does not include any sequencing-depth correction
    

    gX<-g(as.vector(X))
    dim(gX)<-dim(X)
    CVX<-rep(0,dim(X)[2])
    for(j in seq_len(dim(X)[2])){
        CVX[j]=ECVar(X[,j])
    }

    vargx<-stats::var(gX)
    
    varf<-stats::var(gX)-diag(CVX)

    return(varf)
}



LinearCorrectedVariance <- function(X){
    v<-dim(X)[1]
    w<-dim(X)[2]

    #centralise X
    Xbar<-colMeans(X)
    Xtilde<-X-rep(1,v)%*%t(Xbar)

    correctvar<-t(Xtilde)%*%Xtilde/(v-1)-diag(Xbar)
    return(correctvar)
}


LinearCorrectedVarianceSeqDepth <- function(X){
    n<-dim(X)[1]
    p<-dim(X)[2]

    d<-rowSums(X)
    Xstd<-X/(d%*%t(rep(1,p)))

    #centralise X
    Xbar<-colMeans(Xstd)
    SXtilde<-Xstd-rep(1,n)%*%t(Xbar)

    correctvar<-stats::cov(Xstd) -diag(colMeans(Xstd/(d%*%t(rep(1,p))))) 
    return(correctvar)
}
