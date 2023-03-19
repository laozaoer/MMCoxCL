Data_object=function(t,X,Z,b,r,beta,gamma,theta,U,H){
  Ht=H(t)
  if(r>0){
    S=(1+r*Ht*exp(sum(beta*X)+sum(gamma*Z)+theta*b))^(-1/r)
  }
  else{
    S=exp(-Ht*exp(sum(beta*X)+sum(gamma*Z)+theta*b))
  }
  
  return(S-U)
}
Generate_singleT=function(X,Z,b,r,beta,gamma,theta,H){
  U=runif(1,0,1)
  tresult=nleqslv::nleqslv(x=0.1,fn=Data_object,X=X,Z=Z,b=b,r=r,beta=beta,gamma=gamma,theta=theta,U=U,H=H)
  return(tresult$x)
}
Generate_T=function(X,Z,b,r,beta,gamma,theta,n,ni,H){
  result=list()
  length(result)=n
  for (i in 1:n) {
    result[[i]]=rep(0,ni[i])
    for (j in 1:ni[i]) {
      result[[i]][j]=Generate_singleT(X[[i]][j,],Z[i,],b[i],r,beta,gamma,theta,H)
    }
  }
  return(result)
}
data_for_est=function(r,beta,gamma,theta,n,H){
  betadim=length(beta)
  gammadim=length(gamma)
  Z=matrix(runif(n*gammadim,-1,1),nrow = n,ncol = gammadim)
  mi=rep(0,n)
  b=rnorm(n,0,1)
  for(i in 1:n){
    mi[i]=extraDistr::rtpois(1,exp(1.7),a=1,b=8)
  }
  C=list()
  length(C)=n+1
  for(i in 1:n){
    C[[i]]=runif(mi[i],0,1)
  }
  X=list()
  length(X)=n
  for (i in 1:n) {
    X[[i]]=matrix(runif(mi[i]*betadim,-1,1),nrow = mi[i],ncol=betadim)
  }
  
  
  
  rawC=suppressWarnings(Generate_T(X,Z,b,r,beta,gamma,theta,n,mi,H))
  lowC=0
  upC=quantile(as.numeric(as.character(unlist(rawC))),probs = 0.85)
  
  for(i in 1:n){
    C[[i]]=runif(mi[i],lowC,upC)
  }
  C[[n+1]]=upC
  Delta=list()
  length(Delta)=n
  for (i in 1:n) {
    Delta[[i]]=rep(0,mi[i])
    for (j in 1:mi[i]) {
      if(rawC[[i]][j]<=C[[i]][j]){
        Delta[[i]][j]=1
      }
    }
  }
  
  
  return(list(X=X,Z=Z,n=n,ni=mi,r=r,Delta=Delta,C=C))
}
SIM_DATA=function(beta,theta,n,H,k){
  r=0
  gamma=c(0)
  betadim=length(beta)
  gammadim=length(gamma)
  Z=matrix(runif(n*gammadim,-1,1),nrow = n,ncol = gammadim)
  mi=rep(0,n)
  b=rnorm(n,0,1)
  for(i in 1:n){
    mi[i]=sample(c(1:k),1)
  }
  C=list()
  length(C)=n
  for(i in 1:n){
    C[[i]]=runif(mi[i],0,1)
  }
  X=list()
  length(X)=n
  for (i in 1:n) {
    X[[i]]=cbind(runif(mi[i],-1,1),rbinom(mi[i],1,0.5))
  }
  
  
  
  rawC=suppressWarnings(Generate_T(X,Z,b,r,beta,gamma,theta,n,mi,H))
  lowC=quantile(as.numeric(as.character(unlist(rawC))),probs = 0.45)
  upC=quantile(as.numeric(as.character(unlist(rawC))),probs = 0.9)
  
  
  
  for(i in 1:n){
    C[[i]]=cbind(runif(mi[i],0,lowC),runif(mi[i],lowC,upC))
  }
  Delta=list()
  length(Delta)=n
  for (i in 1:n) {
    Delta[[i]]=matrix(0,nrow=mi[i],2)
    for (j in 1:mi[i]) {
      if(rawC[[i]][j]<=C[[i]][j,1]){
        Delta[[i]][j,1]=1
      }else if(rawC[[i]][j]>C[[i]][j,2]){
        Delta[[i]][j,2]=1
      }
    }
  }
  for (i in 1:n) {
    for(j in 1:mi[i]){
      if(Delta[[i]][j,1]==1){
        C[[i]][j,]=c(0,C[[i]][j,1])
      }else if(Delta[[i]][j,2]==1){
        C[[i]][j,]=c(C[[i]][j,2],Inf)
      }
    }
  }
  data=list(X=X,Z=Z,n=n,ni=mi,r=r,Delta=Delta,C=C)
  result=c()
  for(i in 1:data$n){
    resultI=cbind(rep(i,data$ni[i]),data$C[[i]],data$Delta[[i]],data$X[[i]])
    result=rbind(result,resultI)
  }
  return(result)
}



Cox_est=function(Data,Cluster.ind=TRUE,CIlevel=0.95,Tol=1e-5){
  betadim=ncol(Data)-5
  raw_num=nrow(Data)
  cluster_num=length(unique(Data[,1]))
  NAindex=rowSums(is.na(Data))==0
  if(sum(NAindex)==raw_num){
    printresult=paste('The data frame has', raw_num, 'many observations, with', cluster_num, 'level two units (clusters), and', betadim, 'covariates, and there is no missing value in the data')
  }else{
    diff=raw_num-sum(NAindex)
    Data=Data[NAindex,]
    raw_num=nrow(Data)
    cluster_num=length(unique(Data[,1]))
    printresult=paste('The data frame has', raw_num, 'many observations, with', cluster_num, 'level two units (clusters), and', betadim, 'covariates after removing', diff ,'rows with missing values')
    
  }
  print(printresult)
  
  
  
  
  for(i in 6:ncol(Data)){
    if(is.character(Data[,i])){
      Data[,i]=model.matrix(~Data[,i])[, -1] 
    }
  }
  
  cluster_num=length(unique(Data[,1]))
  if(cluster_num==nrow(Data)){
    Cluster.ind=FALSE
    printresult=paste('Since the number of level one units is one within each level two unit, this is not a clustered data. So all observations are treated as independent and identically distributed')
    print(printresult)
  }
  if(Cluster.ind){
    myD=as.matrix(Data)
    betadim=ncol(myD)-5

    myrules=hermite.h.quadrature.rules(20,normalized=FALSE)
    myrules=as.matrix(myrules[[20]])
    bDep=MainFuncClosd(myD,myrules,Tol)
    TheConst=bDep[[3]]-qchisq(CIlevel,1)/2
    CIest=ComputeCI(myD,myrules,0.001,as.matrix(bDep[[2]][1:betadim]),bDep[[2]][betadim+1],as.matrix(bDep[[1]]),0.001,TheConst)
    result=cbind(bDep[[2]],CIest)
    result=as.data.frame(result)
    col_names=colnames(myD)
    rownames(result)[1:betadim]=col_names[6:ncol(myD)]
    # for(i in 1:betadim){
    #   rownames(result)[i]=paste("beta",i,sep = "_")
    # }
    rownames(result)[betadim+1]="theta"
    colnames(result)=c("Est","LB","UB")
  }
  else{
    myD=as.matrix(Data)
    betadim=ncol(myD)-5
    myrulesInd=matrix(c(0,0.5,0,0.5),2,2,byrow = TRUE)
    bInd=MainFuncClosdInd(myD,myrulesInd,Tol)
    TheConst=bInd[[3]]-qchisq(CIlevel,1)/2
    CIest=ComputeCIInd(myD,myrulesInd,0.001,as.matrix(bInd[[2]][1:betadim]),bInd[[2]][betadim],as.matrix(bInd[[1]]),0.001,TheConst)
    result=cbind(bInd[[2]],CIest)
    result=as.data.frame(result)
    col_names=colnames(myD)
    rownames(result)[1:betadim]=col_names[6:ncol(myD)]
    colnames(result)=c("Est","LB","UB")
  }
  return(result)
}

