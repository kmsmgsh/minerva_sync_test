snowapr = function(cls,x,y,k,reverse=F,dyn=F,chunksize=1)
{
  require(parallel)
  p=ncol(x)
  allcombs=genallcombs(p,k)
  ncombs=length(allcombs)
  clusterExport(cls,"do1pset")
  tasks=if(!reverse)
    seq(1,ncombs,chunksize)else
      seq(ncombs,1,-chunksize)
  if (!dyn){
    out=clusterApply(cls,tasks,dochunk,x,y,allcombs,chunksize)
  } else {
    out=clusterApply(cls,tasks,dochunk,x,y,allcombs,chunksize)
  }
  #a=Reduce(rbind,out)
  Reduce(rbind,out)
  #print("result is")
  #print (a)
}

genallcombs=function(p,k){
  allcombs=list()
  for(i in 1:k){
    tmp=combn(1:p,i)
    allcombs=c(allcombs,matrixtolist(tmp,rc=2))
  }
  allcombs
}
#column rc=2,row rc=1
matrixtolist=function(rc,m){
  if(rc==1){
    Map(function(rownum) m[rownum,],1:nrow(m))
  }else Map(function(colnum) m[,colnum],1:ncol(m))
}

dochunk=function(psetsstart,x,y,allcombs,chunksize){
  ncombs=length(allcombs)
  lasttask=min(psetsstart+chunksize-1,ncombs)
  t(sapply(allcombs[psetsstart:lasttask],do1pset,x,y))
}

do1pset=function(onepset,x,y){
  slm=summary(lm(y~x[,onepset]))     
  n0s=ncol(x)-length(onepset)      #onepset 是一个 list，包含了做回归的变量名
  c(slm$adj.r.squared,onepset,rep(0,n0s)) # 返回值第一个元素是r^2的值，后面是使用的变量，然后用0补齐位置
}

snowtest=function(cls,n,p,k,chunksize=1,dyn=F,rvrs=F){
  gendata(n,p)
  snowapr(cls,x,y,k,rvrs,dyn,chunksize)
}


gendata=function(n,p){
  x<<-matrix(rnorm(n*p),ncol=p)
  y<<- x%*% c(rep(0.5,p))+rnorm(n)
}
initmc=function(nworkers){
  makeCluster(nworkers)
}

library(parallel)
#snowtest(initmc(4),100,4,2)
n=10000;p=20;k=3
system.time(snowtest(c8,n,p,k))
system.time(snowtest(initmc(2),n,p,k))
system.time(snowtest(initmc(4),n,p,k))



