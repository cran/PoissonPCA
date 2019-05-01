

simulateXfromLambda<-function(Lambda){
    n<-dim(Lambda)[1]
    p<-dim(Lambda)[2]
    X<-rep(0,n*p)
    dim(X)<-c(n,p)
    for(i in 1:n){
        while(sum(X[i,])==0){ #ensure all row sums at least 1
            for(j in 1:p){
                if(Lambda[i,j]<1e9){#This should probably vary between machines
                    X[i,j]<-stats::rpois(1,Lambda[i,j])
                }else{
                    X[i,j]<-round(stats::rnorm(1,Lambda[i,j],sqrt(Lambda[i,j])))
                }
            }
        }
    }
    return(X)
}


LinearPoissonRank2SimulateLambda<-function(n,p,muu,muv,muw,sigu,sigv,sigw){
    #Simulates a rank-2 n*p matrix of latent Poisson means.
    #The program simulates a vector v of length n from a normal distribution
    #with mean muv and variance sigv^2, and vectors u and w of length p 
    #from normal distributions with means muu and muw and variances sigu^2 and
    #sigw^2. It then simulates X as a Poisson with mean Lambda=1u^T+vw^T.

    u<-stats::rnorm(p)*sigu+muu
    v<-stats::rnorm(n)*sigv+muv
    w<-stats::rnorm(p)*sigw+muw

#We need w to be orthogonal to 1

    w<-w-sum(w)/p
    
    Lambda<-rep(1,n)%*%t(u)+v%*%t(w)
    return(Lambda)
}


rcompositional<-function(p){
    #Returns random p-dimensional compositional data
    #normally distributed with mean 1/p and variance 1/(3p)
    mu.u<-1/p
    sigma.u<-1/p/3
    u<-stats::rnorm(p,mu.u,sigma.u)
    u<-u+(1-sum(u))/length(u)
    while(min(u)<0){
      u[u<0]<-1/p/4
      u<-u+(1-sum(u))/length(u)
    }  
    return(u)
}
    

LinearPoissonRank2SimulateLambdaCompositional<-function(n,p,w,sigma.v,S.mu,S.sigma){
###simulates compostitional data lambda_c=1u^T+vw^T#####
    #u will be simulated as normal with mean 1/p and standard deviation 1/(3p)
    #then shifted to be compositional
    #w will be simulated as normal with mean mu.w and standard deviation sigma.w
    #v will be simulated  in an interval such that 1u^T+vw^T is non-negative
    #and following a scaled beta distribution such that it has mean 0 and 
    #standard deviation sigma.v
    
    #Simulate compositional u
    u<-rcompositional(p)

  
  #v is a scaled and shifted beta distribution with mean 0 and standard 
  #deviation sigma.v
  #The bounds v.left and v.right are chosen so that 1u^T+vw^T is non-negative.
    v<-rep(0,n)
    v.left<--min((u/w)[u/w>=0])
    v.right<--max((u/w)[u/w<=0])
    mu.beta=-v.left/(v.right-v.left)
    sigma.beta=sigma.v/(v.right-v.left)
    c<-(1+1/mu.beta)
    
    alpha<-(c+sigma.beta^2*(1+c)^2)/(c*sigma.beta^2*(1+c)^2)
    beta<-1+c/(c+sigma.beta^2*(1+c)^2)
    for(v.j in 1:n){
        v[v.j]<-(v.right-v.left)*rbeta(1,alpha,beta)+v.left
    }
  
  #S is sequencing depth. Should be set to at least 5
    S.vector<-rnorm(n,S.mu,S.sigma)
    S.vector[S.vector<=5]<-5
    S<-diag(S.vector)

    lambda<-S%*%(rep(1,times=n)%*%t(u)+v%*%t(w))
  
    lambda[lambda<0]<-0
    return(lambda)
}




