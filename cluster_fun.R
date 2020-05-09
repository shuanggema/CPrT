###can use any coef_gen file###

cluster_fun=function(data,pos)
{
library(MASS)
library(fda)
library(stats)
library(ggplot2)
source("cluster_FDA_functions.R")


###  functional data ###


int1=pos[2]-pos[1]
n=nrow(data)



### seperating 4th of the observed points in to test and train data


stest=round(length(pos)/4)
test=seq(1,length(pos),by=4)
postest=pos[test]

 
y=data[,-test]
ytest=data[,test]
 
  
K=K_det(y,pos,int1)
Bsp=create.bspline.basis(rangeval=c(0,2),norder = 3,nbasis = K)
B1=eval.basis(pos, Bsp, returnMatrix=T)
B=B1[-test,]
Btest=B1[test,]


##least squares soln###

bls=matrix(0,nrow=nrow(y),ncol=ncol(B))
for( pop in 1: n)
{
  bls[pop,] =lm(y[pop,]~0+B)$coefficients
}
 

#### penalized clustering###

lambda=c(0.05,0.1,1)#c(0.0001,0.001,0.01,.1,1,2,3)
rho=1
name=c("lam","obj_pen","M_pen")
res=matrix(0,nrow=length(lambda),ncol=length(name))
colnames(res)=name

for(tun in 1:length(lambda))
{
    
    lam=lambda[tun]
     b=bls
    out=objmin(lam,b,B,rho,y,1:n,adap=1)
    best=out[[1]]
     penres2=cut_er_single_meth2(best,bls,ytest,Btest)
    mest2=length(penres2[[2]])
    obj_pen=penres2[1]
    res[tun,]=unlist(c(lam,obj_pen,mest2))
    #cat("sn:",res[tun,],"\n")
    #write.table(res[tun,],file='out/dum1.csv',row.names=F,append=T,quote=F,sep=',')
    
    
  }
  
  res=data.frame(res)
  lam_ob_pen=lambda[which(res$obj_pen==min(res$obj_pen))[1]]
  out=objmin( lam_ob_pen,b,B,rho,y,1:n,adap=1)
  best=out[[1]]
  penres2=cut_er_single_meth2(best,bls,ytest,Btest)
  clusters=penres2[[2]]
  
### reshaping data for plots 
data_long=dat_reshape(n,clusters,data,pos)
data_long$clusters=as.character(data_long$clusters)
  
  
p1 <- ggplot(data_long, aes(x=pos, y=functions, colour=clusters, group=id)) +
  geom_line() 
#plot(p1)
 
return(list(clusters,p1))


} 
  