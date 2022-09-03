Cox_est=function(Data,Cluster.ind=TRUE,CIlevel=0.95,Tol=1e-5){
  if(Cluster.ind){
    myD=Data
    betadim=ncol(myD)-5
    myrules=hermite.h.quadrature.rules(20,normalized=FALSE)
    myrules=as.matrix(myrules[[20]])
    bDep=MainFuncClosd(myD,myrules,Tol)
    TheConst=bDep[[3]]-qchisq(CIlevel,1)/2
    CIest=ComputeCI(myD,myrules,0.001,as.matrix(bDep[[2]][1:betadim]),bDep[[2]][betadim+1],as.matrix(bDep[[1]]),0.001,TheConst)
    result=cbind(bDep[[2]],CIest)
    result=as.data.frame(result)
    for(i in 1:betadim){
      rownames(result)[i]=paste("beta",i,sep = "_")
    }
    rownames(result)[betadim+1]="theta"
    colnames(result)=c("Est","LB","UB")
  }
  else{
    myD=Data
    betadim=ncol(myD)-5
    myrulesInd=matrix(c(0,0.5,0,0.5),2,2,byrow = TRUE)
    bInd=MainFuncClosdInd(myD,myrulesInd,Tol)
    TheConst=bInd[[3]]-qchisq(CIlevel,1)/2
    CIest=ComputeCIInd(myD,myrulesInd,0.001,as.matrix(bInd[[2]][1:betadim]),bInd[[2]][betadim],as.matrix(bInd[[1]]),0.001,TheConst)
    result=cbind(bInd[[2]],CIest)
    result=as.data.frame(result)
    for(i in 1:betadim){
      rownames(result)[i]=paste("beta",i,sep = "_")
    }
    colnames(result)=c("Est","LB","UB")
  }
  return(result)
}