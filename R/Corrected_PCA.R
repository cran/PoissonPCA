Poisson_Corrected_PCA<-function(X,k=dim(X)[2]-1,transformation=NULL,seqdepth=FALSE){
    call=match.call()
    X<-as.matrix(X)
    n<-dim(X)[1]
    if(is.null(transformation)|identical(transformation,"linear")){
        if(seqdepth){
            Sigma<-LinearCorrectedVarianceSeqDepth(X)
        }else{
            Sigma<-LinearCorrectedVariance(X)
        }
        eig<-eigen(Sigma)
        Pcomps<-eig$vectors[,seq_len(k)]
        Pscores<-X%*%Pcomps
        ans<-list("sdev"=sqrt(eig$values),"loadings"=Pcomps,"center"=colMeans(X),"scale"=rep(1,n),"n.obs"=n,"scores"=Pscores,"call"=call,"variance"=Sigma)
        class(ans)<-"princomp"
        return(ans)
    }else{
        if(identical(transformation,"log")){
            transformation<-makelogtransformation(3,6)
            ##a=3 N=6 seems a reasonable default for log transformation
        }
        if(!is.null(transformation$ECVar)){
            Sigma<-TransformedVarianceECV(X,transformation$g,transformation$ECVar)
        }else{
            Sigma<-TransformedVariance(X,transformation$g,transformation$CVar)
        }
        if(any(is.infinite(Sigma))|any(is.nan(Sigma))){
            print("Before Compositional Conversion, Sigma has infinite or NaN values.")
            return(Sigma)
        }

        SigmaStart<-Sigma
        if(seqdepth=="compositional"){
            Sigma<-make_compositional_variance(Sigma)
        }else if(seqdepth=="minvar"){
            Sigma<-make_compositional_min_var(Sigma)
        }else if(seqdepth!=FALSE){
            warning("Sequencing depth method not recognised. No sequencing depth performed.")
            seqdepth=FALSE
        }
        if(any(is.infinite(Sigma))|any(is.nan(Sigma))){
            base::warning("After Compositional Conversion, Sigma has infinite or NaN values.")
            return(Sigma)
        }
        eig<-eigen(Sigma)
        Pscores<-matrix(0,dim(X)[1],k)
        Pmeans<-X
        gX<-transformation$g(as.vector(X))
        dim(gX)<-dim(X)
        mu<-colMeans(gX)        
        V<-eig$vectors
        p<-dim(X)[2]
        if(seqdepth!=FALSE){
            if(k==p-1){
                k=p-2
            }
            proj<-t(V)%*%rep(1/sqrt(p),p)
            index<-which(abs(proj)==max(abs(proj)))
            V<- V - rep(1/sqrt(p),p)%*%t(proj)
            V<- cbind(rep(1/sqrt(p),p),V[,-index])
            for(i in seq_len(p-1)[-1]){
                V[,i]<-V[,i]/sqrt(sum(V[,i]^2))
                V[,(i+1):p]<-V[,(i+1):p]-V[,i]%*%t(V[,i])%*%V[,(i+1):p]
            }
            V[,p]<-V[,p]/sqrt(sum(V[,p]^2))        
            ##Now V should be orthogonal
            D<-c(1,eig$values[-index])
            if(transformation$type=="log"){
                for(i in seq_len(dim(X)[1])){
                    gsl<-get_scores_log(X[i,],V,D,k+1,mu)
                    Pscores[i,]<-gsl$scores[seq_len(k)+1]
                    Pmeans[i,]<-gsl$means
                }
            }else{
                for(i in seq_len(dim(X)[1])){
                    gs<-get_scores(X[i,],V,D,k+1,transformation,mu)
                    Pscores[i,]<-gs$scores[seq_len(k)+1]
                    Pmeans[i,]<-gs$means
                }
            }
        }else{
            D<-c(eig$values)
            if(transformation$type=="log"){
                for(i in seq_len(dim(X)[1])){
                    gsl<-get_scores_log(X[i,],V,D,k,mu)
                    Pscores[i,]<-gsl$scores
                    Pmeans[i,]<-gsl$means
                }
            }else{
                for(i in seq_len(dim(X)[1])){
                    gs<-get_scores(X[i,],V,D,k,transformation,mu)
                    Pscores[i,]<-gs$scores
                    Pmeans[i,]<-gs$means
                }
            }            
        }
        ans<-list("sdev"=eig$values[seq_len(k)],"loadings"=eig$vectors[,seq_len(k)],"center"=mu,"scale"=rep(1,n),"n.obs"=n,"scores"=Pscores,"means"=Pmeans,"variance"=Sigma,"non_compositional_variance"=SigmaStart,"call"=call)
        class(ans)<-c("princomp","transformedprincomp")
        return(ans)
    }
}
