
matrixapply<-function(lambda,func){
    n<-dim(lambda)[1]
    p<-dim(lambda)[2]
    answer<-rep(0,n*p)
    dim(answer)<-c(n,p)
    for(i in 1:n){
        for(j in 1:p){
            answer[i,j]<-func(lambda[i,j])
        }
    }
    return(answer)        
}

setupSimulationAnswer<-function(ColumnNames,sizev,simulation.time,p){
    numsampsize<-length(sizev)
    random.seed<-rep(0,simulation.time*numsampsize)
    dim(random.seed) <- c(numsampsize,simulation.time)

    NumMethods<-length(ColumnNames)
    result_full_mean<-matrix(0,nrow = numsampsize,ncol = NumMethods)
    result_full_sd<-matrix(0,nrow = numsampsize,ncol = NumMethods)
    colnames(result_full_mean)<-ColumnNames
    colnames(result_full_sd)<-ColumnNames
    rownames(result_full_mean)<-rep(0,numsampsize)
    rownames(result_full_sd)<-rep(0,numsampsize)

    for (i in 1:numsampsize){
        for(j in 1:simulation.time){
            random.seed[i,j]<-runif(1, min = 0, max = 10^5)
        }
        n<-sizev[i]
        RowName<-paste0("result",n,"_",p)
        rownames(result_full_mean)[i]<-RowName
        rownames(result_full_sd)[i]<-RowName
    }
    return(list("mean"=result_full_mean,"sd"=result_full_sd,"seed"=random.seed))
}


runSimulations<-function(mu.w,sigma.w,sigma.v,S.sigma,S.mu,sizev,p,transformation,simulation.time=100){
#Simulates transformed lambda as rank-2 1u^T+vw^T where u and w are 
#normally distributed with u being compositional    
#applies random sequencing depth S to the data, normally distributed. 
#p is the dimension of the data    
#transformation is a closure with 4 functions: f, g, solve and CVar
#f evaluates the transformation
#solve reverses the transformation
#g computes an estimator for f(lambda)
#CVar computes an estimator for the conditional variance of g

    ColumnNames<-c("transformed","corrected","compositionaltransformed","Lambdacompositional")

    numsampsize<-length(sizev)
    NumMethods<-length(ColumnNames)

    answer<-setupSimulationAnswer(ColumnNames,sizev,simulation.time,p)
    

    for (i in 1:numsampsize){
        n<-sizev[i]
        result.onesize<-matrix(c(0:0),nrow=simulation.time,ncol=NumMethods)
        
        for(j in 1:simulation.time){
        #save seed for each simulation
            set.seed(answer$random.seed[i,j])
        #w should be centred but not scaled.
            w<-stats::rnorm(p,mu.w,sigma.w)
            w<-scale(w,center = T,scale = F)

            lambda_transform<-LinearPoissonRank2SimulateLambdaCompositional(n,p,w,sigma.v,S.mu,S.sigma)

            lambda<-matrixapply(lambda_transform,transformation$solve)

            X<-simulateXfromLambda(lambda)            
            

            Xtransform<-matrixapply(X,transformation$f)

            
            d<-rowSums(X)
            Xcompositional<-X/(d%*%t(rep(1,p)))
            Xcompositionaltransform<-matrixapply(Xcompositional,transformation$f)

            
            w<-w/((sum(w^2))^0.5)

            result.onesize[j,1]<-abs(sum(eigen(cov(Xtransform))$vectors[,1]*w))
            Sigma<-TransformedVariance(X,transformation$g,transformation$CVar)

            
            result.onesize[j,2]<-abs(sum(eigen(Sigma)$vectors[,1]*w))
            result.onesize[j,3]<-abs(sum(eigen(cov(Xcompositionaltransform))$vectors[,1]*w))
            result.onesize[j,4]<-abs(sum(eigen(cov(lambda_transform))$vectors[,1]*w))
        }
        answer$mean[i,]<-colMeans(result.onesize)
        newcol<-rep(0,NumMethods)
        for(j in 1:NumMethods){
            newcol[j]<-sd(result.onesize[,j])/sqrt(simulation.time)
        }
        answer$sd[i,]<-newcol
    }
    return(answer)
}


runSimulationslinear<-function(mu.w,sigma.w,sigma.v,S.sigma,S.mu,sizev,p,simulation.time=100){

    ColumnNames<-c("originalPCA","linear","linearSeq","Xcompositional","Lambdacompositional")

    numsampsize<-length(sizev)
    NumMethods<-length(ColumnNames)

    answer<-setupSimulationAnswer(ColumnNames,numsampsize,simulation.time)
    

    for (i in 1:numsampsize){
        n<-sizev[i]
        
        result.onesize<-matrix(c(0:0),nrow=simulation.time,ncol=NumMethods)
    
        for(j in 1:simulation.time){
        #save seed for each simulation
            set.seed(answer$random.seed[i,j])
        
    #w should be centred but not scaled.
            w<-stats::rnorm(p,mu.w,sigma.w)
            w<-scale(w,center = T,scale = F)

            lambda<-LinearPoissonRank2SimulateLambdaCompositional(n,p,w,sigma.v,S.mu,S.sigma)
      
        
            X<-simulateXfromLambda(lambda)

            lambdacompositional<-lambda/(rowSums(lambda)%*%t(rep(1,p)))


            d<-rowSums(X)
            Xcompositional<-X/(d%*%t(rep(1,p)))

            w<-w/((sum(w^2))^0.5)

            
            result.onesize[j,1]<-abs(sum(eigen(cov(X))$vectors[,1]*w))
        
            Sigma<-LinearCorrectedVariance(X)

            result.onesize[j,2]<-abs(sum(eigen(Sigma)$vectors[,1]*w))

            Sigma<-LinearCorrectedVarianceSeqDepth(X)

            
            result.onesize[j,3]<-abs(sum(eigen(Sigma)$vectors[,1]*w))
            
            result.onesize[j,4]<-abs(sum(eigen(cov(Xcompositional))$vectors[,1]*w))

            
            result.onesize[j,5]<-abs(sum(eigen(cov(lambdacompositional))$vectors[,1]*w))
        }
    
        answer$mean[i,]<-colMeans(result.onesize)
        newcol<-rep(0,NumMethods)
        for(j in 1:NumMethods){
            newcol[j]<-sd(result.onesize[,j])/sqrt(simulation.time)
        }
        answer$sd[i,]<-newcol
    }

    return(answer)
}


#
#
#mu.w<-0
#sigma.w<-1
#
#sigma.v<-1
#S.sigma<-250
#S.mu<-500
#
#
#sizev<-c(100,500,1000,2000)
#simulation.time<-100
#p<-10

#runSimulations(mu.w,sigma.w,sigma.v,S.sigma,S.mu,sizev,p)


