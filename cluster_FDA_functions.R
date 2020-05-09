
K_det=function(y,pos1,int1)
{
  tmp=1:ncol(y)
  ntest=round(ncol(y)/5)
  test=sample(1:ncol(y),ntest)
  train=tmp[!tmp %in% test]
  ytest=y[,test]
  ytrain=y[,train]
  postest=pos1[test]
  postrain=pos1[train]
  K=(4:10)
  K_opt=NULL
  for(i in 1: nrow(y))
  {
    pred_er=NULL
    for( tun in 1: length(K))
    {
      Bsp=create.bspline.basis(rangeval=c(0,2),norder = 3,nbasis =as.numeric(K[tun]))
      B1=eval.basis(pos1, Bsp, returnMatrix=T)
      Btr=B1[train,]
      Btest=B1[test,]
      bls=lm(ytrain[i,]~0+Btr)$coefficients
      # bls_int=bls[1]
      ypred=apply(Btest,1,function(v)crossprod(bls,v))
      pred_er[tun]=mean((ytest-ypred)^2)
      
    }
    K_opt[i]=K[which(pred_er==min(pred_er,na.rm=T))]
    
  }
  
  K=max(K_opt)
  
  return(K)
}






check=function(v)
{
  if(v[1]>=v[2]) {return(0)}
  if(v[1]<v[2]) {return(1)}
}



b_ADMM=function(ygp,B,i,v,rho,z,bgp,ngp)
{
  tmp1=t(B)%*%B+(ngp-1)*rho*diag(1,ncol(B))
  sm=which(v$Var1==i)
  lr=which(v$Var2==i)
  if(length(sm)>0)
  {
    vtmp1=-apply(v[sm,1:ncol(B)],2,sum)
    ztmp1=apply(z[sm,1:ncol(B)],2,sum)
  }
  if(length(sm)==0)
  {
    ztmp1=vtmp1=rep(0,ncol(B))
  }
  if(length(lr)>0)
  {
    vtmp2=apply(v[lr,1:ncol(B)],2,sum)
    ztmp2=-apply(z[lr,1:ncol(B)],2,sum)
  }
  if(length(lr)==0){ztmp2=vtmp2=rep(0,ncol(B))}
  tmp2=t(B)%*%ygp[i,]+vtmp1+vtmp2+rho*(apply(matrix(bgp[-i,],ncol=ncol(B)),2,sum))+rho*(ztmp1+ztmp2)
  return(solve(tmp1)%*%tmp2)
  
}

z_ADMM=function(lam,rho,v,b,bls,adap)
{
  
  ztmp1=matrix(0,nrow=nrow(v),ncol=ncol(b))
  ztmp1=cbind(ztmp1,v[,(ncol(b)+1):(ncol(b)+2)])
  for( d in 1:nrow(v))
  {
    i1=ztmp1$Var1[d]
    j1=ztmp1$Var2[d]
    if(adap==1)
    {
      lamij= lam/(sum((bls[i1,]-bls[j1,])^2))^(1/2)
    }
    if(adap==2)
    {
      lamij= lam/(sum((b[i1,]-b[j1,])^2))^(1/2)
    }
    if(adap==0)
    {
      lamij= lam
    }
    r=b[i1,]-b[j1,]+v[d,1:ncol(b)]/rho
    norm=(sum(r^2))^(1/2)
    ztmp1[d,1:ncol(b)]=max(0,norm-(lamij/rho))*r/norm
    
  }
  
  return(ztmp1)
}

v_ADMM=function(v,rho,b,z)
{
  vtmp1=matrix(0,nrow=nrow(z),ncol=ncol(b))
  vtmp1=cbind(vtmp1,z[,(ncol(b)+1):(ncol(b)+2)])
  for( d in 1:nrow(z))
  {
    i1=vtmp1$Var1[d]
    j1=vtmp1$Var2[d]
    vtmp1[d,1:ncol(b)]=v[d,1:ncol(b)]+rho*(b[i1,]-b[j1,]-z[d,1:ncol(b)])
  }
  return(vtmp1)
}



objmin=function(lam,blsgp,B,rho,ygp,indx,adap)
{
  ngp=length(indx)
  track1=1:ngp
  track2=1:ngp
  track=expand.grid(track1,track2)
  dum1=apply(track,1,function(v) check(v))
  track=track[which(dum1==1),]
  z=matrix(0,nrow=nrow(track),ncol=ncol(B))
  z=v=data.frame(z,track)
  b_up=t(sapply(1:ngp,function(i) b_ADMM(ygp,B,i,v,rho,z,blsgp,ngp)))
  z_up=z_ADMM(lam,rho,v,b_up,blsgp,adap)
  v_up=v_ADMM(v,rho,b_up,z_up)
  dist=NULL
  for( i in 1: nrow(blsgp))
  {
    dum1=(blsgp[i,]-b_up[i,])^2
    dist[i]=sqrt(sum(dum1))
  }
  iter=1
  diff=NULL
  diff[iter]=max(dist)
  
  while(diff[iter]>0.051)
  {
    b=b_up
    z=z_up
    v=v_up
    b_up=t(sapply(1:ngp,function(i) b_ADMM(ygp,B,i,v,rho,z,b,ngp)))
    z_up=z_ADMM(lam,rho,v,b_up,blsgp,adap)
    v_up=v_ADMM(v,rho,b_up,z_up)
    dist=double(nrow(blsgp))
    for( i in 1: nrow(blsgp))
    {
      dum1=(b[i,]-b_up[i,])^2
      dist[i]=sqrt(sum(dum1))
    }
    iter=iter+1
    diff[iter]=max(dist)
    #print(diff[iter])
    
  }
  result=list(b_up,z_up)
  return(result)
}



betadiff=function(best,indx)
{
  ngp=nrow(best)
  track=expand.grid(1:ngp,1:ngp)
  dum2=apply(track,1,function(v) check(v))
  track=track[which(dum2==1),]
  bdiff=data.frame(diff=NA,track)
  for(d in 1:nrow(track))
  {
    i1=track$Var1[d]
    j1=track$Var2[d]
    dumup=(best[i1,]-best[j1,])^2
    bdiff[d,1]=sqrt(sum(dumup))
  }
  track=expand.grid(indx,indx)
  dum2=apply(track,1,function(v) check(v))
  track=track[which(dum2==1),]
  bdiff[,2:3]=track
  return(bdiff)
}

grouping=function(bdiff,ngp,indx,cutoff)
{
  #cutoff=0.01
  grp=list()
  gps=which(bdiff[,1]<cutoff)
  if(length(gps)>0)
  {
    bdiffgp=matrix(bdiff[gps,],ncol=3)
    ind=1
    d1=1
    indx1=list()
    indx1[[ind]]=d1:length(gps)
    cont="Y"
    
    while(cont=="Y")
    {
      grp[[ind]]=bdiffgp[d1,2:3]
      indx2=NULL
      for( d2 in indx1[[ind]])
      {
        samp2=bdiffgp[d2,2:3]
        if(length(intersect(grp[[ind]],samp2))>0)
        {grp[[ind]]=c(unlist(grp[[ind]]),samp2)}
        if(length(intersect(grp[[ind]],samp2))==0)
        {indx2=c(indx2,d2)}
        
        if(length(indx2)>0)
        {
          for( i in 1:length(indx2))
          {
            if(length(intersect(grp[[ind]],bdiffgp[indx2[i],2:3]))>0)
            {
              if(length(indx2==1)){indx2=NULL}
              if(length(indx2>1)){indx2=indx2[-i]}
              grp[[ind]]=c(unlist(grp[[ind]]),bdiffgp[indx2[i],2:3])
            }
          }
          
        }
      }
      grp[[ind]]=unique(grp[[ind]])  
      
      if(length(indx2)>0)
      {
        ind=ind+1
        indx1[[ind]]=indx2
        d1=indx1[[ind]][1]
        cont="Y"
      }
      if(length(indx2)==0)
      {
        ind=ind+1
        cont="N"
      }
      
    }
    gropuedele=NULL
    for( l in 1: length(grp))
    {
      gropuedele=c(gropuedele,unlist(grp[[l]]))
      
    }
    if(length(gropuedele)< ngp)
    {
      mat=match(indx,gropuedele)
      for( pop in 1:length(indx[is.na(mat)]))
      {
        grp[[ind]]=indx[is.na(mat)][pop]
        ind=ind+1
      }
    }
  }
  if(length(gps)==0)
  {
    grp=list()
    for( pop in 1:length(indx))
    {
      grp[[pop]]=indx[pop]
    }
  }
  return(grp)
}


obj=function(ytest,Btest,best,bls,group_1)
{
  n=nrow(ytest)
  predls= predup=double(n)
  for( pop in 1:n)
  {
    dum2=dum3=double(nrow(Btest))
    for( k1 in 1:nrow(Btest))
    {
      dum2[k1]=(Btest[k1,]%*%bls[pop,]-ytest[pop,k1])^2
      dum3[k1]=(Btest[k1,]%*%best[pop,]-ytest[pop,k1])^2
    }
    
    predls[pop]=mean(dum2)
    predup[pop]=mean(dum3)
  }
  ngp=length(group_1)
  diffls=diffup=double(ngp)
  for( g in 1: ngp)
  {
    blsgp=bls[group_1[[g]],]
    bestgp=best[group_1[[g]],]
    nele=length(group_1[[g]])
    if(nele>1)
    {
      dum1=0
      track=expand.grid(1:nele,1:nele)
      dum2=apply(track,1,function(v) check(v))
      track=track[which(dum2==1),]
      distls=distup=double(nrow(track))
      for(d in 1:nrow(track))
      {
        i1=track$Var1[d]
        j1=track$Var2[d]
        dumls=(blsgp[i1,]-blsgp[j1,])^2
        dumup=(bestgp[i1,]-bestgp[j1,])^2
        distls[d]=sqrt(sum(dumls))
        distup[d]=sqrt(sum(dumup))
      }
    }
    if(nele==1)
    {distls=sum(blsgp^2)
    distup=sum( bestgp^2)
    }
    
    diffls[g]=mean(distls)
    diffup[g]=mean(distup)
  }
  
  objls=mean(predls)+mean(diffls)#.5*predls+lam*diffls
  objup=mean(predup)+mean(diffup)#.5*predup+lam*diffup
  return(objup)
  
}


cut_er_single_meth2=function(best,bls,ytest,Btest)
{
  n=nrow(best)
  bdiff=as.matrix(betadiff(best,1:n))
  cutoff=seq(0.001,2,by=0.005)
  obj_cut=rep(-1,length(cutoff))
  for(cut in 1:length(cutoff))
  {
    group_1=grouping(bdiff,n,1:n,cutoff=cutoff[cut])
    obj_cut[cut]=obj(ytest,Btest,best,bls,group_1)
  }
  
  cutoff_fin=cutoff[which(obj_cut==min(obj_cut))][1]
  grpfin=grouping(bdiff,n,1:n,cutoff=cutoff_fin)
  objgp=obj(ytest,Btest,best,bls,grpfin)
  grp1=NULL
  nG=length(grpfin)
  grp1=NULL
  for( i in 1:n)
  { for( gp in 1: nG)
  {
    dum1=unlist(grpfin[[gp]])
    if(length(which(dum1==i)>0)){grp1[i]=gp}
  } }
  
  return(list(objgp,grpfin))
}



dat_reshape=function(n,clusters,data,pos)
{
  group=rep(0,n)
  for(i in 1: length(clusters))
  {
    grpi=clusters[[i]]
    group[grpi]=i
  }
  
  data=cbind(data,group)
  id=fun=grp=pos1=NULL
  for( i in 1: nrow(data))
  {
    id=c(id,rep(i,length(pos)))
    pos1=c(pos,pos1)
    fun=c(fun,data[i,-ncol(data)])
    grp=c(grp,rep(data[i,ncol(data)],length(pos)))
    
  }
  data_long=cbind(id,pos,fun,grp)
  data_long=data.frame(id=id,pos=pos,functions=fun,clusters=grp)  
  return(data_long) 
}
  
  
  
  






