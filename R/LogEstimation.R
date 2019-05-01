
makelogtransformation<-function(a,N,uselog=6,unbiassed=TRUE){
    force(a)
    force(N)
    force(uselog)
    rangeupper<-200
    gvals0<-log(a)-sum(1/(1:N))
    gvals<-rep(0,rangeupper)
    for(i in seq_len(uselog)){
        gvals[i]<-gvals0
        term<--N*i/a
        for(j in 1:N){
            gvals[i]<-gvals[i]-term
            term<-term*(-1)*(i-j)*(N-j)/a*j/(j+1)^2
        }
    }
    ref<--1
    mode(ref)<-"integer"
    ref<-.C('precalculate',ref,as.double(a),as.integer(N),as.integer(uselog+1))

    
    hh<-NULL
    for(i in uselog:rangeupper){
        gvals[i]=log(i)
    }
    f<-function(x){
        return(log(x))
    }
    solvefunction<-function(x){
        return(exp(x))
    }

    preparehh<-function(mxx){
        gvals<-rep(0,mxx+1)
        hh<<-rep(0,mxx+1)
        agvals<-rep(0,mxx+1)
        for(i in 1:(mxx+1)){
            agvals[i]<-g(i-1)
        }
        for(i in 0:mxx){
            hh[i+1]<<-sum((agvals[1:(i+1)]/(gamma(1:(i+1))))*(agvals[(i+1):1]/(gamma(1:(i+1)))))/(agvals[i+1])^2
        }
    }
    g<-function(x){
        n<-length(x)
        answer<-vector("double",n)
        returnval<-.C('log_g',answer,as.integer(n),as.integer(x),as.integer(ref))
        answer<-returnval[[1]]
        return(answer)
    }
    if(unbiassed){
        ECVar<-function(x){
###Needs a vector valued x
            n<-length(x)

            answers<-vector("double",1)
            returnval<-.C('log_ECVar_unbiassed',answers,as.integer(n),as.integer(ref),as.integer(x),as.double(a),as.integer(N),as.integer(uselog+1))
            answers<-returnval[[1]]
            return(answers)
        }
    }else{
        ECVar<-function(x){
###Needs a vector valued x
            n<-length(x)
            answers<-vector("double",1)
            returnval<-.C('log_ECVar',answers,as.integer(n),as.double(a),as.integer(N),as.integer(x))
            answers<-returnval[[1]]
            return(answers)
        }
    }
    answer<-list("type"="log","f"=f,"g"=g,"solve"=solvefunction,"ECVar"=ECVar)
    return(answer)
}



