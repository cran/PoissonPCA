library(PoissonPCA)

n<-1000   #number of samples
p<-10   #dimension of samples
r<-3    #rank of covariance matrix

Z<-stats::rnorm(n*r)
U<-stats::rnorm(p*r)/sqrt(r)
dim(Z)<-c(n,r)
dim(U)<-c(r,p)
Lambda<-exp(Z%*%U)
Sigma<-t(U)%*%U
seqdepth<-exp(stats::rnorm(n)+5) #sequencing depth is log-normal
LambdaObserved<-diag(seqdepth)%*%Lambda
X<-stats::rpois(n*p,as.vector(LambdaObserved))
dim(X)<-c(n,p)


logtransform<-makelogtransformation(3,4)


system.time(correctedPCA<-Poisson_Corrected_PCA(X,2,logtransform,seqdepth="minvar"))

#Compare the components to their projections onto the true space.
(t(solve(U%*%t(U),U%*%correctedPCA$components))%*%U)%*%correctedPCA$components

M<-t(solve(U%*%t(U),U%*%correctedPCA$components))

plot((M%*%U)[1,],correctedPCA$components[,1])
plot((M%*%U)[2,],correctedPCA$components[,2])
#confirm that principal components are MU

#so now log(Lambda)=mu+ZU
#Projection should be mu+ZM^{-1}MU
#so corresponding weights should be ZM^{-1}


truescores<-Z%*%t(solve(M%*%t(M),M))#this should be the correct value of the scores

#Compare the scores to the observed scores
#plot((Z%*%)[,1],correctedPCA$scores[,1])
#abline(0,5)
plot((Z%*%solve(U%*%t(U),U%*%correctedPCA$components))[,2],correctedPCA$scores[,2])
