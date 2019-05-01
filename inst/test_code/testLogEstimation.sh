#!/bin/bash
samplesize=$1
filename=tableofvalues_$1.txt
ggraphname=meang_$1.pdf
CVgraphname=varg_$1.pdf


R --no-save<<EOF
library(PoissonPCA)
uvals<-(1:100)/20-1
sampsize<-$samplesize
logtransform<-makelogtransformation(3,4)

CVarx<-rep(0,1000)
gx<-rep(0,1000)

meang<-rep(0,100)
varg<-rep(0,100)
meanCV<-rep(0,100)

for(j in 1:100){
    u<-uvals[j]
    psample<-rpois(1000,exp(u))
    for(i in 1:1000){
        gx[i]<-logtransform\$g(psample[i])
    }
    meang[j]=mean(gx)
    varg[j]=var(gx)
    meanCV[j]=logtransform\$ECVar(psample)
}
pdf("$ggraphname")
plot(uvals,meang,type='l',xlab="u",ylab="E(g(x))")
abline(0,1)
dev.off()

pdf("$CVgraphname")
plot(uvals,varg,type='l',xlab="u",ylab="Var(g(x))")
points(uvals,meanCV,type='l',col="red")
dev.off()

write.table(cbind(uvals,meang,varg,meanCV),"$filename")
EOF

