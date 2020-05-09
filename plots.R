###can use any coef_gen file###
library(ggplot2)
library(MASS)

source("src/main_1_func.R")
source("src/settings.R")
source("src/coef_gen3.R")
source("src/coef_gen4.R")
source("src/curve_plot_func.R")
source("src/coef_gen5.R")

set1=setting()
randID=runif(1,1,10)

setin=1

delta=set1$delta[setin];error=set1$error[setin];m=set1$m[setin]; shapes=set1[setin,1:m];int1=set1$int1[setin]; sd=set1$sd[setin]; cdist=set1$cdist[setin]; coefgen=set1$coefgen[setin];  adap=set1$adap[setin];  ni=set1[setin,(m+1):(2*m)]
error=0.2
ni=as.matrix(ni,ncol=1)
  n=sum(ni)
  pos1=seq(0,2,by=int1)[-1]
  y1=matrix(0,nrow=n,ncol=(length(pos1)+1))
  y1_noerror=matrix(0,nrow=m,ncol=length(pos1))
  
  ind=0
  
  for(gp in 1:m)
  {
      y1tmp=ygen4(shapes[gp],ni[gp],error,int1)
      y1_noerror[gp,]=ygen4_noerr_pl(1,shapes[gp],int1)
      #y1tmp_noerror[gp,]=ygen4_noerr_pl(1,shapes[gp],int1)
    
    y1[(ind+1):(ind+ni[gp]),]=cbind(y1tmp,rep(gp,ni[gp]))
   # y1_noerror=rbind(y1_noerror,y1tmp_noerror)
    plot(pos1,apply(y1tmp,2,mean),xlab=gp,type="l")
    ind=ind+ni[gp]
  }
  
  
data=NULL
id=fun=grp=pos=NULL
for( i in 1: nrow(y1))
{
  id=c(id,rep(i,length(pos1)))
  pos=c(pos,pos1)
  fun=c(fun,y1[i,-ncol(y1)])
  grp=c(grp,rep(y1[i,ncol(y1)],length(pos1)))
  
}

data=cbind(id,pos,fun,grp)
data=data.frame(id=id,pos=pos,fun=fun,grp=grp)  
  
  head(data)
  

###Group 1
clustr=1
shapei=y1_noerror[clustr,]
data_sub=data[which(data$grp==clustr),]
data_sub=data_sub[,-4]
data_sub=data.frame(data_sub,smo=0)



tmp1=data.frame(id=rep(0,length(pos1)),pos=pos1,fun=shapei,smo=1)

data_sub=rbind(data_sub,tmp1)
data_sub$id=as.character(data_sub$id)
data_sub$smo=as.character(data_sub$smo)

head(data_sub)


p1 <- ggplot(data_sub, aes(x=pos, y=fun,group=id)) +
ylab(" ")+
geom_line(aes(linetype=id))+
scale_linetype_manual(values=c( "solid","longdash","dotted", "twodash", "dotdash"))+
ggtitle("Cluster 1")+
theme(plot.title = element_text(hjust = 0.5))+
 theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
plot(p1)





###Group 2
clustr=2
shapei=y1_noerror[clustr,]
data_sub=data[which(data$grp==clustr),]
data_sub=data_sub[,-4]
data_sub=data.frame(data_sub,smo=0)



tmp1=data.frame(id=rep(0,length(pos1)),pos=pos1,fun=shapei,smo=1)

data_sub=rbind(data_sub,tmp1)
data_sub$id=as.character(data_sub$id)
data_sub$smo=as.character(data_sub$smo)



p2 <- ggplot(data_sub, aes(x=pos, y=fun,group=id)) +
  ylab(" ")+
  geom_line(aes(linetype=id))+
  scale_linetype_manual(values=c( "solid","longdash","dotted", "twodash", "dotdash"))+
  ggtitle("Cluster 2")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")

plot(p2)





###Group 3
clustr=3
shapei=y1_noerror[clustr,]
data_sub=data[which(data$grp==clustr),]
data_sub=data_sub[,-4]
data_sub=data.frame(data_sub,smo=0)



tmp1=data.frame(id=rep(0,length(pos1)),pos=pos1,fun=shapei,smo=1)

data_sub=rbind(data_sub,tmp1)
data_sub$id=as.character(data_sub$id)
data_sub$smo=as.character(data_sub$smo)



p3 <- ggplot(data_sub, aes(x=pos, y=fun,group=id)) +
  ylab(" ")+
  geom_line(aes(linetype=id))+
  scale_linetype_manual(values=c( "solid","longdash","dotted", "twodash", "dotdash"))+
  ggtitle("Cluster 3")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
plot(p3)



###Group 4
clustr=4
shapei=y1_noerror[clustr,]
data_sub=data[which(data$grp==clustr),]
data_sub=data_sub[,-4]
data_sub=data.frame(data_sub,smo=0)



tmp1=data.frame(id=rep(0,length(pos1)),pos=pos1,fun=shapei,smo=1)

data_sub=rbind(data_sub,tmp1)
data_sub$id=as.character(data_sub$id)
data_sub$smo=as.character(data_sub$smo)



p4 <- ggplot(data_sub, aes(x=pos, y=fun,group=id)) +
  ylab(" ")+
  geom_line(aes(linetype=id))+
  scale_linetype_manual(values=c( "solid","longdash","dotted", "twodash", "dotdash"))+
  ggtitle("Cluster 4")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
plot(p4)




###Group 5
clustr=5
shapei=y1_noerror[clustr,]
data_sub=data[which(data$grp==clustr),]
data_sub=data_sub[,-4]
data_sub=data.frame(data_sub,smo=0)



tmp1=data.frame(id=rep(0,length(pos1)),pos=pos1,fun=shapei,smo=1)

data_sub=rbind(data_sub,tmp1)
data_sub$id=as.character(data_sub$id)
data_sub$smo=as.character(data_sub$smo)



p5 <- ggplot(data_sub, aes(x=pos, y=fun,group=id)) +
  ylab(" ")+
  geom_line(aes(linetype=id))+
  scale_linetype_manual(values=c( "solid","longdash","dotted", "twodash", "dotdash"))+
  ggtitle("Cluster 5")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
plot(p5)


## Group 6

clustr=6
shapei=y1_noerror[clustr,]
data_sub=data[which(data$grp==clustr),]
data_sub=data_sub[,-4]
data_sub=data.frame(data_sub,smo=0)



tmp1=data.frame(id=rep(0,length(pos1)),pos=pos1,fun=shapei,smo=1)

data_sub=rbind(data_sub,tmp1)
data_sub$id=as.character(data_sub$id)
data_sub$smo=as.character(data_sub$smo)



p6 <- ggplot(data_sub, aes(x=pos, y=fun,group=id)) +
  ylab(" ")+
  geom_line(aes(linetype=id))+
  scale_linetype_manual(values=c( "solid","longdash","dotted", "twodash", "dotdash"))+
  ggtitle("Cluster 6")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
plot(p6)





## Group 7

clustr=7
shapei=y1_noerror[clustr,]
data_sub=data[which(data$grp==clustr),]
data_sub=data_sub[,-4]
data_sub=data.frame(data_sub,smo=0)



tmp1=data.frame(id=rep(0,length(pos1)),pos=pos1,fun=shapei,smo=1)

data_sub=rbind(data_sub,tmp1)
data_sub$id=as.character(data_sub$id)
data_sub$smo=as.character(data_sub$smo)



p7 <- ggplot(data_sub, aes(x=pos, y=fun,group=id)) +
  ylab(" ")+
  geom_line(aes(linetype=id))+
  scale_linetype_manual(values=c( "solid","longdash","dotted", "twodash", "dotdash"))+
  ggtitle("Cluster 7")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
plot(p7)




## Group 8

clustr=8
shapei=y1_noerror[clustr,]
data_sub=data[which(data$grp==clustr),]
data_sub=data_sub[,-4]
data_sub=data.frame(data_sub,smo=0)



tmp1=data.frame(id=rep(0,length(pos1)),pos=pos1,fun=shapei,smo=1)

data_sub=rbind(data_sub,tmp1)
data_sub$id=as.character(data_sub$id)
data_sub$smo=as.character(data_sub$smo)



p8 <- ggplot(data_sub, aes(x=pos, y=fun,group=id)) +
  ylab(" ")+
  geom_line(aes(linetype=id))+
  scale_linetype_manual(values=c( "solid","longdash","dotted", "twodash", "dotdash"))+
  ggtitle("Cluster 8")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
plot(p8)


multiplot(p1, p2, p3, p4, p5,p6,p7,p8, cols=4)






