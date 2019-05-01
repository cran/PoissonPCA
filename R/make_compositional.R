
make_compositional_variance<-function(Sigma){
#Converts a variance matrix Sigma to a composition variance matrix
    p<-dim(Sigma)[1] 
    if(dim(Sigma)[2]!=p){#check Sigma is square
        stop("Variance matrix must be square.")
    }
    a<-solve(diag(rep(p,p))+matrix(1,nrow=p,ncol=p),rowSums(Sigma))
    return(Sigma-rep(1,p)%*%t(a)-a%*%t(rep(1,p)))
}

make_compositional_min_var<-function(Sigma){
#subtracts the largest possible constant matrix from Sigma while remaining non-negative definite
    p<-dim(Sigma)[1] 
    if(dim(Sigma)[2]!=p){#check Sigma is square
        stop("Variance matrix must be square.")
    }
    eig<-eigen(Sigma)
    teig<-sum(eig$values[-1][eig$values[-1]>0])
    neig<-sum(cumsum(eig$values[-1][eig$values[-1]>0])<teig*0.9)+2
###take largest eigenvalues covering more than 90% of variance.
    thr<-eig$values[neig]
    Vgood<-eig$vectors[,seq_len(neig)]
    VVo<-rep(1,p)-Vgood%*%(t(Vgood)%*%rep(1,p))
    VVo<-VVo/sqrt(sum(VVo^2))
###Projection of 1 onto space orthogonal to V.
    neweval<-rep(thr,p)
    neweval[seq_len(neig)]<-eig$values[seq_len(neig)]-thr
    Sigmathr<-diag(rep(thr,p))+Vgood%*%diag(neweval[seq_len(neig)])%*%t(Vgood)#-VVo%*%t(VVo)*thr
    Sigmathr<-(Sigmathr+t(Sigmathr))/2  #Sigmathr could be numerically not symmetric
    Sigma_c<-make_compositional_variance(Sigmathr)
    Sigma_star<-Sigma_c+matrix(1/p,nrow=p,ncol=p)
    sigma_ssq<-det(Sigmathr)/p/det(Sigma_star)
    return(Sigmathr-sigma_ssq*matrix(1,nrow=p,ncol=p))
}


make_compositional_min_var_old<-function(Sigma){
#subtracts the largest possible constant matrix from Sigma while remaining non-negative definite
    p<-dim(Sigma)[1] 
    if(dim(Sigma)[2]!=p){#check Sigma is square
        stop("Variance matrix must be square.")
    }
    Sigma_c<-make_compositional_variance(Sigma)
    Sigma_star<-Sigma_c+matrix(1/p,nrow=p,ncol=p)
    sigma_ssq<-det(Sigma)/p/det(Sigma_star)
    return(Sigma-sigma_ssq*matrix(1,nrow=p,ncol=p))
}
