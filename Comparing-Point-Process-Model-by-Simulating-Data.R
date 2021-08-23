library(spatstat)
#######################################################################################
#Pemodelan Simulate Data#
#######################################################################################
comp.model<-function(kappa,r,mu,rep){
    nclust <- function(x0, y0, radius, n) {
      return(runifdisc(n, radius, centre=c(x0, y0)))
    }
   cp<- rPoissonCluster(kappa, 0.2, nclust, radius=r, n=mu, win=square(2),nsim = rep)
   modp<-lapply(cp, ppm)
   mod11<-lapply(cp, function(i){kppm(i~1,"Thomas", method="clik2")})
   mod12<-lapply(cp, function(j){kppm(j~1,"Thomas", method="palm")})
   mod21<-lapply(cp, function(k){kppm(k~1,"MatClust", method="clik2")})
   mod22<-lapply(cp, function(ii){kppm(ii~1,"MatClust", method="palm")})
   mod31<-lapply(cp, function(ij){kppm(ij~1,"Cauchy", method="clik2")})
   mod32<-lapply(cp, function(ik){kppm(ik~1,"Cauchy", method="palm")})
   mod41<-lapply(cp, function(iii){kppm(iii~1,"VarGamma", method="clik2")})
   mod42<-lapply(cp, function(iij){kppm(iij~1,"VarGamma", method="palm")})
   mod51<-lapply(cp, function(iik){kppm(iik~1,"LGCP", method="clik2")})
   mod52<-lapply(cp, function(ijk){kppm(ijk~1,"LGCP", method="palm")})
   AICp<-sapply(modp, AIC)
   AIC11<-sapply(mod11,AIC)
   AIC12<-sapply(mod12,AIC)
   AIC21<-sapply(mod21,AIC)
   AIC22<-sapply(mod22,AIC)
   AIC31<-sapply(mod31,AIC)
   AIC32<-sapply(mod32,AIC)
   AIC41<-sapply(mod41,AIC)
   AIC42<-sapply(mod42,AIC)
   AIC51<-sapply(mod51,AIC)
   AIC52<-sapply(mod52,AIC)
   BICp<-sapply(modp, BIC)
   BIC11<-sapply(mod11,BIC)
   BIC12<-sapply(mod12,BIC)
   BIC21<-sapply(mod21,BIC)
   BIC22<-sapply(mod22,BIC)
   BIC31<-sapply(mod31,BIC)
   BIC32<-sapply(mod32,BIC)
   BIC41<-sapply(mod41,BIC)
   BIC42<-sapply(mod42,BIC)
   BIC51<-sapply(mod51,BIC)
   BIC52<-sapply(mod52,BIC)
   aa<-rbind(AICp,AIC11,AIC12,AIC21,AIC22,AIC31,AIC32,AIC41,AIC42,AIC51,AIC52)
   dd<-rbind(BICp,BIC11,BIC12,BIC21,BIC22,BIC31,BIC32,BIC41,BIC42,BIC51,BIC52)
   bb<-as.data.frame(aa)
   ee<-as.data.frame(dd)
   cc<-as.matrix(bb)
   ff<-as.matrix(ee)
   colnames(cc)<-NULL
   row.names(cc)<-NULL
   colnames(ff)<-NULL
   row.names(ff)<-NULL
   allaic<-cc
   allbic<-ff
   method<-c("PP","ThomasC","ThomasP","MatClustC","MatClustP","CauchyC","CauchyP","VarGammaC","VarGammaP",
             "LGCPC","LGCPP")
   codea=matrix(0,11,rep)
   codeb=matrix(0,11,rep)
   pr.best=matrix(0,11,2)
   for (j in 1:11){
     for (r in 1:rep){
       if (cc[j,r]==min(cc[,r])){
         codea[j,r]=1}
       else{
         codea[j,r]=0}
       if (ff[j,r]==min(ff[,r])){
         codeb[j,r]=1}
       else{
         codeb[j,r]=0}
   pr.best[j,1]=(sum(codea[j,]))/rep
   pr.best[j,2]=(sum(codeb[j,]))/rep
   rownames(pr.best)=method
   colnames(pr.best)=c("AIC","BIC")
       }
   }
   aAICp<-mean(AICp)
   aAIC11<-mean(AIC11)
   aAIC12<-mean(AIC12)
   aAIC21<-mean(AIC21)
   aAIC22<-mean(AIC22)
   aAIC31<-mean(AIC31)
   aAIC32<-mean(AIC32)
   aAIC41<-mean(AIC41)
   aAIC42<-mean(AIC42)
   aAIC51<-mean(AIC51)
   aAIC52<-mean(AIC52)
   aBICp<-mean(BICp)
   aBIC11<-mean(BIC11)
   aBIC12<-mean(BIC12)
   aBIC21<-mean(BIC21)
   aBIC22<-mean(BIC22)
   aBIC31<-mean(BIC31)
   aBIC32<-mean(BIC32)
   aBIC41<-mean(BIC41)
   aBIC42<-mean(BIC42)
   aBIC51<-mean(BIC51)
   aBIC52<-mean(BIC52)
   resultAIC<-c(aAICp,aAIC11,aAIC12,aAIC21,aAIC22,aAIC31,aAIC32,aAIC41,aAIC42,aAIC51,aAIC52)
   resultBIC<-c(aBICp,aBIC11,aBIC12,aBIC21,aBIC22,aBIC31,aBIC32,aBIC41,aBIC42,aBIC51,aBIC52)
   resultname<-cbind(method,resultAIC,resultBIC)
   colnames(resultname)<-c("Method","aAIC","aBIC")
   for (i in 1:11){
     if (resultAIC[i]==min(resultAIC)||resultBIC[i]==min(resultBIC)){
       min<-resultname[i,]
     }
   }
print(list(result=resultname,min=min, All.AIC=allaic, All.BIC=allbic,powertest=pr.best))
}
s1=comp.model(5,0.1,5,100)
