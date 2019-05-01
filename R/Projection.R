
get_scores<-function(X,V,d,k,transformation,mu){
#X is the data
#The estimated latent variance is S=VDV^T where D is a diagonal matrix 
#with entries d. k is the number of principal components to project onto.
#transformation is the transformation function.

#The objective function is 
#  X^Tlog(Lambda)+sum(Lambda)-(f(Lambda)-mu-P)^TS^{-1}(f(Lambda)-mu-P)

# mu is assumed known
# P is the projection of f(Lambda) onto the first k principal components.
# To simplify the algebra, we will let f(Lambda)=mu+Va for some vector a
# This makes the objective function
#  X^Tlog(Lambda)+sum(Lambda)-sum_{i=k+1}^p {a_i}^2/d_i
# Suppose we let g=f^{-1} and g' be the derivative of g. Then the 
# derivative of this objective function is
# sum (X_j (d{Lambda_j}/d{a_i})/Lambda_j+sum(d{Lambda_j}/d{a_i})-

    p<-length(X)
    aold<-rep(0,p)
    Lambda<-X
    Lambda[Lambda<=0]<-1
    a<-t(V)%*%(log(Lambda)-mu)
    while(sum((a-aold)^2)>1e-15){
        aold<-a
        a[seq_len(k)]<-rep(0,k)
        Lambda<-X+V%*%(a/d)
        Lambda[Lambda<=0]<-1e-4
        a<-t(V)%*%(log(Lambda)-mu)
    }
    return(list("scores"=a,"means"=Lambda))
# The derivative of the objective function with respect to Lambda_i is
#  X_i/Lambda_i+1-2(f'(Lambda_i)-mu-P)^TS^{-1}(f(Lambda)-mu-P)    
}




get_scores_log<-function(X,V,d,k,mu){
#X is the data
#The estimated latent variance is S=VDV^T where D is a diagonal matrix 
#with entries d. k is the number of principal components to project onto.
#transformation is the transformation function.

#The objective function is 
#  X^Tlog(Lambda)+sum(Lambda)-(log(Lambda)-mu-P)^TS^{-1}(log(Lambda)-mu-P)

# mu is assumed known
# P is the projection of f(Lambda) onto the first k principal components.
# To simplify the algebra, we will let f(Lambda)=mu+Va for some vector a
# This makes the objective function
#  X^Tlog(Lambda)+sum(Lambda)-sum_{i=k+1}^p {a_i}^2/d_i
# Suppose we let g=f^{-1} and g' be the derivative of g. Then the 
# derivative of this objective function is
# V(X-Lambda)-2D^{-1}a
# where Lambda=e^{mu+Va}    

# We solve this via Newton-Raphson method.


#    print(length(X))
    
    p<-length(X)
    
    dstar<-d
    dstar[dstar<0.01]<-0.01
    dstar[seq_len(k)]<-0
#    print(dstar)
    M<-V%*%diag(dstar)%*%t(V)
    
    answer<-vector("double",k+p+2)
#    print(dim(M))
    returnval<-.C('get_scores_log_cpp',answer,as.integer(p),as.integer(X),as.double(as.vector(V)),as.double(dstar),as.integer(k),as.double(mu),as.double(as.vector(M)))
    answer<-returnval[[1]]
    
    return(list("scores"=answer[seq_len(k)],"means"=answer[seq_len(p)+k],"convergence"=answer[p+k+1],"accuracy"=answer[p+k+2]))
}


