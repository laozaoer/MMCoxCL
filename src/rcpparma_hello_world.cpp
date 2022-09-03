// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]

arma::umat TmatL(const arma::mat&Inspec,const arma::vec&Timepoints){
    int n=Inspec.n_rows;
    int m=Timepoints.n_elem;
    arma::umat result(m,n);
    for(int i=0;i<n;i++){
        result.col(i)=Timepoints<=Inspec(i,0);
    }
    return result;
}
// [[Rcpp::export]]
arma::umat TmatR(const arma::mat&Inspec,const arma::vec&Timepoints){
    int n=Inspec.n_rows;
    int m=Timepoints.n_elem;
    arma::umat result(m,n);
    for(int i=0;i<n;i++){
        result.col(i)=Timepoints<=Inspec(i,1);
    }
    return result;
}

// [[Rcpp::export]]

arma::umat TmatLR(const arma::mat&Inspec,const arma::vec&Timepoints){
    int n=Inspec.n_rows;
    int m=Timepoints.n_elem;
    arma::umat result(m,n);
    for(int i=0;i<n;i++){
        
        result.col(i)=Timepoints<=Inspec(i,1)&&Timepoints>Inspec(i,0);
    }
    for(int i=0;i<n;i++){
        if(Inspec.row(i).has_inf()){
            result.col(i)=Timepoints<std::pow(10,10);
        }
    }
    return result;
}

// [[Rcpp::export]]
arma::mat elementwise_pow(const arma::mat&A,const arma::mat&p){
    int I=A.n_rows;
    int J=A.n_cols;
    arma::mat result=arma::zeros(I,J);
    for(int i=0;i<I;i++){
        for(int j=0;j<J;j++){
            result(i,j)=std::pow(A(i,j),p(i,j));
        }
        
    }
    return result;
}

// [[Rcpp::export]]

arma::mat Likelihood_rules(const arma::mat&rules,const arma::mat&beta,const double&theta,const arma::mat&gamma,const arma::mat&Data,
                           const arma::umat&tL,const arma::umat&tR){
    arma::mat index=sort(unique(Data.col(0)));
    int order=rules.n_rows;
    arma::mat covariate=Data.submat(0,5,Data.n_rows-1,Data.n_cols-1);
    arma::mat betaX=covariate*beta;
    arma::mat ParPart=exp(repmat(betaX,1,order)+repmat(trans(theta*rules.col(0)),Data.n_rows,1));
    
    arma::mat LambdaR=trans(trans(gamma)*tR);
    arma::mat LambdaL=trans(trans(gamma)*tL);
    arma::mat expL=exp(-repmat(LambdaL,1,order)%ParPart);
    arma::mat expR=exp(-repmat(LambdaR,1,order)%ParPart);
    arma::mat Likelihood=elementwise_pow(1-expR,repmat(Data.col(3),1,order))%
        elementwise_pow(expL-expR,1-repmat(Data.col(3),1,order)-repmat(Data.col(4),1,order))%
        elementwise_pow(expL,repmat(Data.col(4),1,order));
    
    // arma::mat logLikelihood=repmat(Data.col(3),1,order)%log(1-expR)+(1-repmat(Data.col(3),1,order)-repmat(Data.col(4),1,order))%log(expL-expR)+
    //                        repmat(Data.col(4),1,order)%log(expL);
    // arma::mat Likelihood=exp(logLikelihood);
    return Likelihood;
}



// [[Rcpp::export]]

arma::mat WeightFunc(const arma::mat&rules,const arma::mat&beta,const double&theta,const arma::mat&gamma,const arma::mat&Data,
                     const arma::umat&tL,const arma::umat&tR){
    arma::vec normvec=arma::normpdf(rules.col(0));
    arma::vec weightvec=rules.col(1)%exp(pow(rules.col(0),2))%normvec;
    int order=rules.n_rows;
    // arma::mat Oriweight=repmat(trans(weightvec),Data.n_rows,1);
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    arma::mat LikelihoodMat=Likelihood_rules(rules,beta,theta,gamma,Data,tL,tR);
    arma::mat weightMat=arma::zeros(n,order);
    for(int i=1;i<(n+1);i++){
        arma::umat rowindices=find(Data.col(0)==i);
        arma::mat Eachcluster=LikelihoodMat.rows(rowindices);
        arma::mat Rowprod = arma::prod(Eachcluster,0);
        arma::mat Finalweight=Rowprod%trans(weightvec)/accu(Rowprod%trans(weightvec));
        weightMat.row(i-1)=Finalweight;
    }
    return weightMat;
}



// [[Rcpp::export]]
double LogLikeliValue(const arma::mat&rules,const arma::mat&beta,const double&theta,const arma::mat&gamma,const arma::mat&Data,
                      const arma::umat&tL,const arma::umat&tR){
    arma::vec normvec=arma::normpdf(rules.col(0));
    arma::vec weightvec=rules.col(1)%exp(pow(rules.col(0),2))%normvec;
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    arma::mat LikelihoodMat=Likelihood_rules(rules,beta,theta,gamma,Data,tL,tR);
    arma::vec resultvec=arma::zeros(n);
    for(int i=1;i<(n+1);i++){
        arma::umat rowindices=find(Data.col(0)==i);
        arma::mat Eachcluster=LikelihoodMat.rows(rowindices);
        arma::mat Rowprod = arma::prod(Eachcluster,0);
        arma::mat Finalweight=Rowprod*(weightvec);
        resultvec(i-1)=Finalweight(0,0);
    }
    double value=sum(log(resultvec));
    return value;
    
}

// [[Rcpp::export]]
double TargetLogLikeliValue(const arma::vec&par,const arma::mat&rules,const arma::mat&Data){
    
    int betadim=Data.n_cols-5;
    arma::mat tktemp=sort(unique(join_cols(Data.col(1),Data.col(2))));
    arma::vec tk=tktemp.col(0);
    tk=tk.subvec(1,tk.n_elem-2);
    tk=tk(find(tk<std::pow(10,10)));
    double maxL=Data.col(1).max();
    arma::mat Rinspection=Data.col(2);
    Rinspection=Rinspection.rows(find(Rinspection<std::pow(10,10)));
    double maxR=Rinspection.max();
    arma::uvec tkfinity=find(tk<=(maxL)&&tk<=(maxR));
    tk=tk(tkfinity);
    
    arma::mat L=Data.col(1);
    arma::mat R=Data.col(2);
    
    arma::umat tL=TmatL(join_rows(Data.col(1),Data.col(2)),tk);
    arma::umat tR=TmatR(join_rows(Data.col(1),Data.col(2)),tk);

    arma::mat beta=par.subvec(0,betadim-1);
    double theta=par(betadim);
    arma::mat gamma=par.subvec(betadim+1,par.n_elem-1);
    
    arma::vec normvec=arma::normpdf(rules.col(0));
    arma::vec weightvec=rules.col(1)%exp(pow(rules.col(0),2))%normvec;
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    arma::mat LikelihoodMat=Likelihood_rules(rules,beta,theta,gamma,Data,tL,tR);
    arma::vec resultvec=arma::zeros(n);
    for(int i=1;i<(n+1);i++){
        arma::umat rowindices=find(Data.col(0)==i);
        arma::mat Eachcluster=LikelihoodMat.rows(rowindices);
        arma::mat Rowprod = arma::prod(Eachcluster,0);
        arma::mat Finalweight=Rowprod*(weightvec);
        resultvec(i-1)=Finalweight(0,0);
    }
    double value=sum(log(resultvec));
    return -value;
    
}


// [[Rcpp::export]]
double LogLikeliValueInd(const arma::mat&rules,const arma::mat&beta,const double&theta,const arma::mat&gamma,const arma::mat&Data,
                         const arma::umat&tL,const arma::umat&tR){
    arma::vec weightvec=rules.col(1);
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    arma::mat LikelihoodMat=Likelihood_rules(rules,beta,theta,gamma,Data,tL,tR);
    arma::vec resultvec=arma::zeros(n);
    for(int i=1;i<(n+1);i++){
        arma::umat rowindices=find(Data.col(0)==i);
        arma::mat Eachcluster=LikelihoodMat.rows(rowindices);
        arma::mat Rowprod = arma::prod(Eachcluster,0);
        arma::mat Finalweight=Rowprod*(weightvec);
        resultvec(i-1)=Finalweight(0,0);
    }
    double value=sum(log(resultvec));
    return value;
    
}

// [[Rcpp::export]]
arma::field<arma::mat> Updateonce(const arma::mat&rules,const arma::mat&beta,const double&theta,const arma::mat&gamma,const arma::mat&Data,
                                  const arma::umat&tL,const arma::umat&tR,const arma::umat&tLR){
    arma::mat WeightMat=WeightFunc(rules,beta,theta,gamma,Data,tL,tR);
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    int order=rules.n_rows;
    int m=gamma.n_rows;
    arma::mat covariate=Data.submat(0,5,Data.n_rows-1,Data.n_cols-1);
    arma::mat betaX=covariate*beta;
    arma::mat ParPart=exp(repmat(betaX,1,order)+repmat(trans(theta*rules.col(0)),Data.n_rows,1));
    
    arma::mat LambdaR=repmat(trans(trans(gamma)*tR),1,order);
    arma::mat LambdaL=repmat(trans(trans(gamma)*tL),1,order);
    arma::mat LambdaLR=repmat(trans(trans(gamma)*tLR),1,order);
    arma::mat expL=LambdaL%ParPart;
    arma::mat expR=LambdaR%ParPart;
    arma::mat expLR=LambdaLR%ParPart;
    
    arma::mat CubegammaR=repmat(trans(pow(gamma,3)),Data.n_rows,1)%trans(tR);
    arma::mat CubegammaLR=repmat(trans(pow(gamma,3)),Data.n_rows,1)%trans(tLR);
    arma::mat InversgammaL=repmat(trans(pow(gamma,-1)),Data.n_rows,1)%trans(tL);
    
    
    arma::mat SharedTerm1=repmat(Data.col(3),1,order)%pow((1+expR),-1);
    arma::mat SharedTerm2=(1-repmat(Data.col(3),1,order)-repmat(Data.col(4),1,order))%pow((1+expLR),-1);
    arma::mat SharedTerm3=(1-repmat(Data.col(3),1,order))%expL;
    
    arma::mat Term1=trans(SharedTerm1%pow((LambdaR),-1));
    arma::mat Term2=trans(SharedTerm2%pow((LambdaLR),-1));
    arma::mat Term3=trans((1-repmat(Data.col(3),1,order))%ParPart);
    // std::cout<<pow((LambdaLR),-1);
    arma::mat Amat(n,m),Bmat(n,m);
    
    arma::mat FirstBigMat=SharedTerm1+SharedTerm2-SharedTerm3;
    
    arma::mat SecondBigMat=FirstBigMat+2*SharedTerm3;
    
    arma::mat FirstDerivMat(n,beta.n_rows+1);
    arma::mat SecondDeriv=arma::zeros(beta.n_rows+1,beta.n_rows+1);
    arma::mat SecondDerivMat=arma::zeros(beta.n_rows+1,beta.n_rows+1);
    for(int i=0;i<n;i++){
        arma::umat rowindices=find(Data.col(0)==(i+1));
        arma::mat AclusterI=Term1.cols(rowindices)*CubegammaR.rows(rowindices)+
            Term2.cols(rowindices)*CubegammaLR.rows(rowindices);
        arma::mat BclusterI=Term3.cols(rowindices)*InversgammaL.rows(rowindices);
        
        Amat.row(i)=WeightMat.row(i)*AclusterI;
        Bmat.row(i)=WeightMat.row(i)*BclusterI;
        // std::cout<<Amat;
        arma::mat FirstBigMatij=trans(FirstBigMat.rows(rowindices));
        arma::mat SecondBigMatij=SecondBigMat.rows(rowindices);
        arma::mat FirstDerivI(order,beta.n_rows+1);
        for(int r=0;r<order;r++){
            
            arma::mat Xij=covariate.rows(rowindices);
            arma::vec onesvec=arma::ones(Xij.n_rows);
            arma::mat Wijr=join_rows(Xij,rules(r,0)*onesvec*theta);
            FirstDerivI.row(r)=FirstBigMatij.row(r)*Wijr;
            arma::mat SecondDerivI=arma::zeros(beta.n_rows+1,beta.n_rows+1);
            
            for(int j=0;j<Xij.n_rows;j++){
                SecondDerivI=SecondDerivI+SecondBigMatij(j,r)*(trans(Wijr.row(j))*Wijr.row(j));
            }
            SecondDeriv=SecondDeriv-2*WeightMat(i,r)*SecondDerivI;
            
            
        }
        FirstDerivMat.row(i)=WeightMat.row(i)*FirstDerivI;
        
        
    }
    arma::mat FirstDeriv=trans(sum(FirstDerivMat,0));
    SecondDeriv(beta.n_rows,beta.n_rows)=SecondDeriv(beta.n_rows,beta.n_rows)+FirstDeriv(beta.n_rows);
    arma::mat updategamma=pow((trans(sum(Amat,0))%pow(trans(sum(Bmat,0)),-1)),0.25);
    // std::cout<<(trans(sum(Amat,0)));
    arma::umat Extremevalue=find(updategamma<std::pow(10,-10));
    updategamma.rows(Extremevalue)=arma::ones(Extremevalue.n_rows)*std::pow(10,-10);
    
    arma::vec lastpar(beta.n_rows+1);
    lastpar.subvec(0,beta.n_rows-1)=beta;
    lastpar(beta.n_rows)=std::log(theta);
    
    arma::vec updatePar=lastpar-solve(SecondDeriv,FirstDeriv);
    updatePar(beta.n_rows)=std::exp(updatePar(beta.n_rows));
    
    
    arma::field<arma::mat> result(2);
    result(0)=updategamma;
    result(1)=updatePar;
    return result;
    
}

// [[Rcpp::export]]
arma::field<arma::mat> UpdateonceInd(const arma::mat&rules,const arma::mat&beta,const double&theta,const arma::mat&gamma,const arma::mat&Data,
                                     const arma::umat&tL,const arma::umat&tR,const arma::umat&tLR){
    // arma::mat WeightMat=WeightFunc(rules,beta,theta,gamma,Data,tL,tR);
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    arma::mat WeightMat=arma::ones(n,2)*0.5;
    int order=rules.n_rows;
    int m=gamma.n_rows;
    arma::mat covariate=Data.submat(0,5,Data.n_rows-1,Data.n_cols-1);
    arma::mat betaX=covariate*beta;
    arma::mat ParPart=exp(repmat(betaX,1,order)+repmat(trans(theta*rules.col(0)),Data.n_rows,1));
    
    arma::mat LambdaR=repmat(trans(trans(gamma)*tR),1,order);
    arma::mat LambdaL=repmat(trans(trans(gamma)*tL),1,order);
    arma::mat LambdaLR=repmat(trans(trans(gamma)*tLR),1,order);
    arma::mat expL=LambdaL%ParPart;
    arma::mat expR=LambdaR%ParPart;
    arma::mat expLR=LambdaLR%ParPart;
    
    arma::mat CubegammaR=repmat(trans(pow(gamma,3)),Data.n_rows,1)%trans(tR);
    arma::mat CubegammaLR=repmat(trans(pow(gamma,3)),Data.n_rows,1)%trans(tLR);
    arma::mat InversgammaL=repmat(trans(pow(gamma,-1)),Data.n_rows,1)%trans(tL);
    
    
    arma::mat SharedTerm1=repmat(Data.col(3),1,order)%pow((1+expR),-1);
    arma::mat SharedTerm2=(1-repmat(Data.col(3),1,order)-repmat(Data.col(4),1,order))%pow((1+expLR),-1);
    arma::mat SharedTerm3=(1-repmat(Data.col(3),1,order))%expL;
    
    arma::mat Term1=trans(SharedTerm1%pow((LambdaR),-1));
    arma::mat Term2=trans(SharedTerm2%pow((LambdaLR),-1));
    arma::mat Term3=trans((1-repmat(Data.col(3),1,order))%ParPart);
    // std::cout<<pow((LambdaLR),-1);
    arma::mat Amat(n,m),Bmat(n,m);
    
    arma::mat FirstBigMat=SharedTerm1+SharedTerm2-SharedTerm3;
    
    arma::mat SecondBigMat=FirstBigMat+2*SharedTerm3;
    
    arma::mat FirstDerivMat(n,beta.n_rows);
    arma::mat SecondDeriv=arma::zeros(beta.n_rows,beta.n_rows);
    arma::mat SecondDerivMat=arma::zeros(beta.n_rows,beta.n_rows);
    for(int i=0;i<n;i++){
        arma::umat rowindices=find(Data.col(0)==(i+1));
        arma::mat AclusterI=Term1.cols(rowindices)*CubegammaR.rows(rowindices)+
            Term2.cols(rowindices)*CubegammaLR.rows(rowindices);
        arma::mat BclusterI=Term3.cols(rowindices)*InversgammaL.rows(rowindices);
        
        Amat.row(i)=WeightMat.row(i)*AclusterI;
        Bmat.row(i)=WeightMat.row(i)*BclusterI;
        // std::cout<<Amat;
        arma::mat FirstBigMatij=trans(FirstBigMat.rows(rowindices));
        arma::mat SecondBigMatij=SecondBigMat.rows(rowindices);
        arma::mat FirstDerivI(order,beta.n_rows);
        for(int r=0;r<order;r++){
            
            arma::mat Xij=covariate.rows(rowindices);
            arma::vec onesvec=arma::ones(Xij.n_rows);
            arma::mat Wijr=Xij;
            FirstDerivI.row(r)=FirstBigMatij.row(r)*Wijr;
            arma::mat SecondDerivI=arma::zeros(beta.n_rows,beta.n_rows);
            
            for(int j=0;j<Xij.n_rows;j++){
                SecondDerivI=SecondDerivI+SecondBigMatij(j,r)*(trans(Wijr.row(j))*Wijr.row(j));
            }
            SecondDeriv=SecondDeriv-2*WeightMat(i,r)*SecondDerivI;
            
            
        }
        FirstDerivMat.row(i)=WeightMat.row(i)*FirstDerivI;
        // std::cout<<WeightMat.row(i);
        
    }
    arma::mat FirstDeriv=trans(sum(FirstDerivMat,0));
    // SecondDeriv(beta.n_rows,beta.n_rows)=SecondDeriv(beta.n_rows,beta.n_rows)+FirstDeriv(beta.n_rows);
    arma::mat updategamma=pow((trans(sum(Amat,0))%pow(trans(sum(Bmat,0)),-1)),0.25);
    // std::cout<<(trans(sum(Amat,0)));
    arma::umat Extremevalue=find(updategamma<std::pow(10,-10));
    updategamma.rows(Extremevalue)=arma::ones(Extremevalue.n_rows)*std::pow(10,-10);
    
    arma::vec lastpar(beta.n_rows);
    lastpar.subvec(0,beta.n_rows-1)=beta;
    
    arma::vec updatePar=lastpar-solve(SecondDeriv,FirstDeriv);
    
    
    arma::field<arma::mat> result(2);
    result(0)=updategamma;
    result(1)=updatePar;
    // std::cout<<result(1);
    return result;
    
}




// [[Rcpp::export]]
arma::field<arma::mat> MainFuncClosd(const arma::mat&Data,const arma::mat&rules,const double&Tol){
    int betadim=Data.n_cols-5;
    arma::mat tktemp=sort(unique(join_cols(Data.col(1),Data.col(2))));
    arma::vec tk=tktemp.col(0);
    tk=tk.subvec(1,tk.n_elem-2);
    // tk=tk(find(tk<std::pow(10,10)));
    
    double maxL=Data.col(1).max();
    arma::mat Rinspection=Data.col(2);
    Rinspection=Rinspection.rows(find(Rinspection<std::pow(10,10)));
    double maxR=Rinspection.max();
    arma::uvec tkfinity=find(tk<=(maxL)&&tk<=(maxR));
    tk=tk(tkfinity);
    // if(any(Data.col(1)==tk.max())){
    //     tk=tk.subvec(0,tk.n_elem-2);
    // }
    
    
    arma::mat L=Data.col(1);
    arma::mat R=Data.col(2);
    
    arma::umat tL=TmatL(join_rows(Data.col(1),Data.col(2)),tk);
    arma::umat tR=TmatR(join_rows(Data.col(1),Data.col(2)),tk);
    arma::umat tLR=TmatLR(join_rows(Data.col(1),Data.col(2)),tk);
    arma::mat beta0=arma::ones(betadim,1)*0.5;
    arma::mat gamma0=arma::ones(tk.n_elem,1)*0.1;
    double theta0=0.5;
    double absdiff=1000;
    arma::mat betatheta;
    // arma::mat gammaold=gamma0;
    // arma::mat betaold=beta0;
    // double thetaold=theta0;
    // arma::mat Totalold=join_cols(gammaold,betaold,arma::ones(1,1)*thetaold);
    arma::field<arma::mat> newresult;
    int iter=1;
    do{ 
        arma::mat oldbeta=beta0;
        double oldtheta=theta0;
        arma::mat oldgamma=gamma0;
        newresult=Updateonce(rules,beta0,theta0,gamma0,Data,tL,tR,tLR);
        gamma0=newresult(0);
        betatheta=newresult(1);
        beta0=betatheta.submat(0,0,betadim-1,0);
        theta0=betatheta(betadim,0);
        arma::mat Totalnew=newresult(1);
        arma::mat Totalold=join_cols(oldbeta,arma::ones(1,1)*oldtheta);
        absdiff=accu(abs((Totalold-Totalnew)%pow(Totalold,-1)));
        iter=iter+1;
        // std::cout<<iter;
        // std::cout<<absdiff;
        
    } while (iter<500&absdiff>Tol);
    double LogLike=LogLikeliValue(rules,beta0,theta0,gamma0,Data,tL,tR);
    
    arma::field<arma::mat> finalresult(4);
    finalresult(0)=newresult(0);
    finalresult(1)=newresult(1);
    finalresult(2)=arma::ones(1,1)*LogLike;
    finalresult(3)=arma::ones(1,1)*iter;
    return finalresult;
}

// [[Rcpp::export]]
arma::field<arma::mat> MainFuncClosdInd(const arma::mat&Data,const arma::mat&rules,const double&Tol){
    int betadim=Data.n_cols-5;
    arma::mat tktemp=sort(unique(join_cols(Data.col(1),Data.col(2))));
    arma::vec tk=tktemp.col(0);
    tk=tk.subvec(1,tk.n_elem-2);
    // tk=tk(find(tk<std::pow(10,10)));
    
    double maxL=Data.col(1).max();
    arma::mat Rinspection=Data.col(2);
    Rinspection=Rinspection.rows(find(Rinspection<std::pow(10,10)));
    double maxR=Rinspection.max();
    arma::uvec tkfinity=find(tk<=(maxL)&&tk<=(maxR));
    tk=tk(tkfinity);

    
    
    arma::mat L=Data.col(1);
    arma::mat R=Data.col(2);
    
    arma::umat tL=TmatL(join_rows(Data.col(1),Data.col(2)),tk);
    arma::umat tR=TmatR(join_rows(Data.col(1),Data.col(2)),tk);
    arma::umat tLR=TmatLR(join_rows(Data.col(1),Data.col(2)),tk);
    arma::mat beta0=arma::ones(betadim,1)*0.5;
    arma::mat gamma0=arma::ones(tk.n_elem,1)*0.1;
    double theta0=0.5;
    double absdiff=1000;
    arma::mat betatheta;

    arma::field<arma::mat> newresult;
    int iter=1;
    do{ 
        arma::mat oldbeta=beta0;
        arma::mat oldgamma=gamma0;
        newresult=UpdateonceInd(rules,beta0,theta0,gamma0,Data,tL,tR,tLR);
        gamma0=newresult(0);
        betatheta=newresult(1);
        beta0=betatheta.submat(0,0,betadim-1,0);
        arma::mat Totalnew=newresult(1);
        arma::mat Totalold=oldbeta;
        absdiff=accu(abs((Totalold-Totalnew)%pow(Totalold,-1)));
        iter=iter+1;

        
    } while (iter<500&absdiff>Tol);
    double LogLike=LogLikeliValueInd(rules,beta0,theta0,gamma0,Data,tL,tR);
    
    arma::field<arma::mat> finalresult(4);
    finalresult(0)=newresult(0);
    finalresult(1)=newresult(1);
    finalresult(2)=arma::ones(1,1)*LogLike;
    finalresult(3)=arma::ones(1,1)*iter;
    return finalresult;
}

// [[Rcpp::export]]
arma::field<arma::mat> UpdateonceProfileBeta(const arma::uword&Indicator,const arma::mat&rules,const arma::mat&beta,const double&theta,const arma::mat&gamma,const arma::mat&Data,
                                             const arma::umat&tL,const arma::umat&tR,const arma::umat&tLR){
    arma::mat WeightMat=WeightFunc(rules,beta,theta,gamma,Data,tL,tR);
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    int order=rules.n_rows;
    int m=gamma.n_rows;
    arma::mat covariate=Data.submat(0,5,Data.n_rows-1,Data.n_cols-1);
    
    arma::mat Deletecovariate=covariate;
    Deletecovariate.shed_col(Indicator-1);
    arma::mat betaX=covariate*beta;
    arma::mat ParPart=exp(repmat(betaX,1,order)+repmat(trans(theta*rules.col(0)),Data.n_rows,1));
    arma::mat LambdaR=repmat(trans(trans(gamma)*tR),1,order);
    arma::mat LambdaL=repmat(trans(trans(gamma)*tL),1,order);
    arma::mat LambdaLR=repmat(trans(trans(gamma)*tLR),1,order);
    arma::mat expL=LambdaL%ParPart;
    arma::mat expR=LambdaR%ParPart;
    arma::mat expLR=LambdaLR%ParPart;
    
    arma::mat CubegammaR=repmat(trans(pow(gamma,3)),Data.n_rows,1)%trans(tR);
    arma::mat CubegammaLR=repmat(trans(pow(gamma,3)),Data.n_rows,1)%trans(tLR);
    arma::mat InversgammaL=repmat(trans(pow(gamma,-1)),Data.n_rows,1)%trans(tL);
    
    
    arma::mat SharedTerm1=repmat(Data.col(3),1,order)%pow((1+expR),-1);
    arma::mat SharedTerm2=(1-repmat(Data.col(3),1,order)-repmat(Data.col(4),1,order))%pow((1+expLR),-1);
    arma::mat SharedTerm3=(1-repmat(Data.col(3),1,order))%expL;
    
    arma::mat Term1=trans(SharedTerm1%pow((LambdaR),-1));
    arma::mat Term2=trans(SharedTerm2%pow((LambdaLR),-1));
    arma::mat Term3=trans((1-repmat(Data.col(3),1,order))%ParPart);
    // std::cout<<pow((LambdaLR),-1);
    arma::mat Amat(n,m),Bmat(n,m);
    
    arma::mat FirstBigMat=SharedTerm1+SharedTerm2-SharedTerm3;
    
    arma::mat SecondBigMat=FirstBigMat+2*SharedTerm3;
    
    arma::mat FirstDerivMat(n,beta.n_rows);
    arma::mat SecondDeriv=arma::zeros(beta.n_rows,beta.n_rows);
    arma::mat SecondDerivMat=arma::zeros(beta.n_rows,beta.n_rows);
    
    for(int i=0;i<n;i++){
        arma::umat rowindices=find(Data.col(0)==(i+1));
        arma::mat AclusterI=Term1.cols(rowindices)*CubegammaR.rows(rowindices)+
            Term2.cols(rowindices)*CubegammaLR.rows(rowindices);
        arma::mat BclusterI=Term3.cols(rowindices)*InversgammaL.rows(rowindices);
        
        Amat.row(i)=WeightMat.row(i)*AclusterI;
        Bmat.row(i)=WeightMat.row(i)*BclusterI;
        // std::cout<<Amat;
        arma::mat FirstBigMatij=trans(FirstBigMat.rows(rowindices));
        arma::mat SecondBigMatij=SecondBigMat.rows(rowindices);
        arma::mat FirstDerivI(order,beta.n_rows);
        
        for(int r=0;r<order;r++){
            
            arma::mat Xij=Deletecovariate.rows(rowindices);
            arma::vec onesvec=arma::ones(Xij.n_rows);
            arma::mat Wijr=join_rows(Xij,rules(r,0)*onesvec*theta);
            FirstDerivI.row(r)=FirstBigMatij.row(r)*Wijr;
            arma::mat SecondDerivI=arma::zeros(beta.n_rows,beta.n_rows);
            
            for(int j=0;j<Xij.n_rows;j++){
                SecondDerivI=SecondDerivI+SecondBigMatij(j,r)*(trans(Wijr.row(j))*Wijr.row(j));
            }
            SecondDeriv=SecondDeriv-2*WeightMat(i,r)*SecondDerivI;
            
            
        }
        FirstDerivMat.row(i)=WeightMat.row(i)*FirstDerivI;
        
        
    }
    arma::mat FirstDeriv=trans(sum(FirstDerivMat,0));
    SecondDeriv(beta.n_rows-1,beta.n_rows-1)=SecondDeriv(beta.n_rows-1,beta.n_rows-1)+FirstDeriv(beta.n_rows-1);
    arma::mat updategamma=pow((trans(sum(Amat,0))%pow(trans(sum(Bmat,0)),-1)),0.25);
    arma::umat Extremevalue=find(updategamma<std::pow(10,-10));
    updategamma.rows(Extremevalue)=arma::ones(Extremevalue.n_rows)*std::pow(10,-10);
    
    arma::vec lastpar(beta.n_rows);
    arma::mat Deletebeta=beta;
    Deletebeta.shed_row(Indicator-1);
    arma::vec updatePar;
    if(beta.n_rows==1){
        lastpar(beta.n_rows-1)=std::log(theta);
        updatePar=lastpar-solve(SecondDeriv,FirstDeriv);
        updatePar(beta.n_rows-1)=std::exp(updatePar(beta.n_rows-1));
        updatePar.insert_rows(Indicator-1,beta.row(Indicator-1));
    }else{
        
        
        lastpar.subvec(0,beta.n_rows-2)=Deletebeta;
        lastpar(beta.n_rows-1)=std::log(theta);
        
        updatePar=lastpar-solve(SecondDeriv,FirstDeriv);
        updatePar(beta.n_rows-1)=std::exp(updatePar(beta.n_rows-1));
        
        updatePar.insert_rows(Indicator-1,beta.row(Indicator-1));
    }
    arma::field<arma::mat> result(2);
    result(0)=updategamma;
    result(1)=updatePar;
    return result;
    
}


// [[Rcpp::export]]
arma::field<arma::mat> UpdateonceProfileBetaInd(const arma::uword&Indicator,const arma::mat&rules,const arma::mat&beta,const double&theta,const arma::mat&gamma,const arma::mat&Data,
                                             const arma::umat&tL,const arma::umat&tR,const arma::umat&tLR){
    arma::mat WeightMat=WeightFunc(rules,beta,theta,gamma,Data,tL,tR);
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    int order=rules.n_rows;
    int m=gamma.n_rows;
    arma::mat covariate=Data.submat(0,5,Data.n_rows-1,Data.n_cols-1);
    
    arma::mat Deletecovariate=covariate;
    Deletecovariate.shed_col(Indicator-1);
    arma::mat betaX=covariate*beta;
    arma::mat ParPart=exp(repmat(betaX,1,order)+repmat(trans(theta*rules.col(0)),Data.n_rows,1));
    arma::mat LambdaR=repmat(trans(trans(gamma)*tR),1,order);
    arma::mat LambdaL=repmat(trans(trans(gamma)*tL),1,order);
    arma::mat LambdaLR=repmat(trans(trans(gamma)*tLR),1,order);
    arma::mat expL=LambdaL%ParPart;
    arma::mat expR=LambdaR%ParPart;
    arma::mat expLR=LambdaLR%ParPart;
    
    arma::mat CubegammaR=repmat(trans(pow(gamma,3)),Data.n_rows,1)%trans(tR);
    arma::mat CubegammaLR=repmat(trans(pow(gamma,3)),Data.n_rows,1)%trans(tLR);
    arma::mat InversgammaL=repmat(trans(pow(gamma,-1)),Data.n_rows,1)%trans(tL);
    
    
    arma::mat SharedTerm1=repmat(Data.col(3),1,order)%pow((1+expR),-1);
    arma::mat SharedTerm2=(1-repmat(Data.col(3),1,order)-repmat(Data.col(4),1,order))%pow((1+expLR),-1);
    arma::mat SharedTerm3=(1-repmat(Data.col(3),1,order))%expL;
    
    arma::mat Term1=trans(SharedTerm1%pow((LambdaR),-1));
    arma::mat Term2=trans(SharedTerm2%pow((LambdaLR),-1));
    arma::mat Term3=trans((1-repmat(Data.col(3),1,order))%ParPart);
    // std::cout<<pow((LambdaLR),-1);
    arma::mat Amat(n,m),Bmat(n,m);
    
    arma::mat FirstBigMat=SharedTerm1+SharedTerm2-SharedTerm3;
    
    arma::mat SecondBigMat=FirstBigMat+2*SharedTerm3;
    
    arma::mat FirstDerivMat(n,beta.n_rows-1);
    arma::mat SecondDeriv=arma::zeros(beta.n_rows-1,beta.n_rows-1);
    arma::mat SecondDerivMat=arma::zeros(beta.n_rows-1,beta.n_rows-1);
    
    for(int i=0;i<n;i++){
        arma::umat rowindices=find(Data.col(0)==(i+1));
        arma::mat AclusterI=Term1.cols(rowindices)*CubegammaR.rows(rowindices)+
            Term2.cols(rowindices)*CubegammaLR.rows(rowindices);
        arma::mat BclusterI=Term3.cols(rowindices)*InversgammaL.rows(rowindices);
        
        Amat.row(i)=WeightMat.row(i)*AclusterI;
        Bmat.row(i)=WeightMat.row(i)*BclusterI;
        // std::cout<<Amat;
        arma::mat FirstBigMatij=trans(FirstBigMat.rows(rowindices));
        arma::mat SecondBigMatij=SecondBigMat.rows(rowindices);
        arma::mat FirstDerivI(order,beta.n_rows-1);
        
        for(int r=0;r<order;r++){
            
            arma::mat Xij=Deletecovariate.rows(rowindices);
            arma::vec onesvec=arma::ones(Xij.n_rows);
            arma::mat Wijr=Xij;
            FirstDerivI.row(r)=FirstBigMatij.row(r)*Wijr;
            arma::mat SecondDerivI=arma::zeros(beta.n_rows-1,beta.n_rows-1);
            
            for(int j=0;j<Xij.n_rows;j++){
                SecondDerivI=SecondDerivI+SecondBigMatij(j,r)*(trans(Wijr.row(j))*Wijr.row(j));
            }
            SecondDeriv=SecondDeriv-2*WeightMat(i,r)*SecondDerivI;
            
            
        }
        FirstDerivMat.row(i)=WeightMat.row(i)*FirstDerivI;
        
        
    }
    arma::mat FirstDeriv=trans(sum(FirstDerivMat,0));
    // SecondDeriv(beta.n_rows-1,beta.n_rows-1)=SecondDeriv(beta.n_rows-1,beta.n_rows-1)+FirstDeriv(beta.n_rows-1);
    arma::mat updategamma=pow((trans(sum(Amat,0))%pow(trans(sum(Bmat,0)),-1)),0.25);
    arma::umat Extremevalue=find(updategamma<std::pow(10,-10));
    updategamma.rows(Extremevalue)=arma::ones(Extremevalue.n_rows)*std::pow(10,-10);
    
    arma::vec lastpar(beta.n_rows-1);
    arma::mat Deletebeta=beta;
    Deletebeta.shed_row(Indicator-1);
    arma::vec updatePar;

    if(beta.n_rows==1){

        updatePar=beta;
    }else{
        
        
        lastpar.subvec(0,beta.n_rows-2)=Deletebeta;

        updatePar=lastpar-solve(SecondDeriv,FirstDeriv);

        updatePar.insert_rows(Indicator-1,beta.row(Indicator-1));
    }
    arma::field<arma::mat> result(2);
    result(0)=updategamma;
    result(1)=updatePar;
    return result;
    
}



// [[Rcpp::export]]
arma::field<arma::mat> UpdateonceProfileTheta(const arma::mat&rules,const arma::mat&beta,const double&theta,const arma::mat&gamma,const arma::mat&Data,
                                              const arma::umat&tL,const arma::umat&tR,const arma::umat&tLR){
    arma::mat WeightMat=WeightFunc(rules,beta,theta,gamma,Data,tL,tR);
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    int order=rules.n_rows;
    int m=gamma.n_rows;
    arma::mat covariate=Data.submat(0,5,Data.n_rows-1,Data.n_cols-1);
    arma::mat betaX=covariate*beta;
    arma::mat ParPart=exp(repmat(betaX,1,order)+repmat(trans(theta*rules.col(0)),Data.n_rows,1));
    
    arma::mat LambdaR=repmat(trans(trans(gamma)*tR),1,order);
    arma::mat LambdaL=repmat(trans(trans(gamma)*tL),1,order);
    arma::mat LambdaLR=repmat(trans(trans(gamma)*tLR),1,order);
    arma::mat expL=LambdaL%ParPart;
    arma::mat expR=LambdaR%ParPart;
    arma::mat expLR=LambdaLR%ParPart;
    
    arma::mat CubegammaR=repmat(trans(pow(gamma,3)),Data.n_rows,1)%trans(tR);
    arma::mat CubegammaLR=repmat(trans(pow(gamma,3)),Data.n_rows,1)%trans(tLR);
    arma::mat InversgammaL=repmat(trans(pow(gamma,-1)),Data.n_rows,1)%trans(tL);
    
    
    arma::mat SharedTerm1=repmat(Data.col(3),1,order)%pow((1+expR),-1);
    arma::mat SharedTerm2=(1-repmat(Data.col(3),1,order)-repmat(Data.col(4),1,order))%pow((1+expLR),-1);
    arma::mat SharedTerm3=(1-repmat(Data.col(3),1,order))%expL;
    
    arma::mat Term1=trans(SharedTerm1%pow((LambdaR),-1));
    arma::mat Term2=trans(SharedTerm2%pow((LambdaLR),-1));
    arma::mat Term3=trans((1-repmat(Data.col(3),1,order))%ParPart);
    // std::cout<<pow((LambdaLR),-1);
    arma::mat Amat(n,m),Bmat(n,m);
    
    arma::mat FirstBigMat=SharedTerm1+SharedTerm2-SharedTerm3;
    
    arma::mat SecondBigMat=FirstBigMat+2*SharedTerm3;
    
    arma::mat FirstDerivMat(n,beta.n_rows);
    arma::mat SecondDeriv=arma::zeros(beta.n_rows,beta.n_rows);
    arma::mat SecondDerivMat=arma::zeros(beta.n_rows,beta.n_rows);
    for(int i=0;i<n;i++){
        arma::umat rowindices=find(Data.col(0)==(i+1));
        arma::mat AclusterI=Term1.cols(rowindices)*CubegammaR.rows(rowindices)+
            Term2.cols(rowindices)*CubegammaLR.rows(rowindices);
        arma::mat BclusterI=Term3.cols(rowindices)*InversgammaL.rows(rowindices);
        
        Amat.row(i)=WeightMat.row(i)*AclusterI;
        Bmat.row(i)=WeightMat.row(i)*BclusterI;
        // std::cout<<Amat;
        arma::mat FirstBigMatij=trans(FirstBigMat.rows(rowindices));
        arma::mat SecondBigMatij=SecondBigMat.rows(rowindices);
        arma::mat FirstDerivI(order,beta.n_rows);
        for(int r=0;r<order;r++){
            
            arma::mat Xij=covariate.rows(rowindices);
            arma::vec onesvec=arma::ones(Xij.n_rows);
            arma::mat Wijr=Xij;
            FirstDerivI.row(r)=FirstBigMatij.row(r)*Wijr;
            arma::mat SecondDerivI=arma::zeros(beta.n_rows,beta.n_rows);
            
            for(int j=0;j<Xij.n_rows;j++){
                SecondDerivI=SecondDerivI+SecondBigMatij(j,r)*(trans(Wijr.row(j))*Wijr.row(j));
            }
            SecondDeriv=SecondDeriv-2*WeightMat(i,r)*SecondDerivI;
            
            
        }
        FirstDerivMat.row(i)=WeightMat.row(i)*FirstDerivI;
        
        
    }
    arma::mat FirstDeriv=trans(sum(FirstDerivMat,0));
    // SecondDeriv(beta.n_rows,beta.n_rows)=SecondDeriv(beta.n_rows,beta.n_rows)+FirstDeriv(beta.n_rows);
    arma::mat updategamma=pow((trans(sum(Amat,0))%pow(trans(sum(Bmat,0)),-1)),0.25);
    arma::umat Extremevalue=find(updategamma<std::pow(10,-10));
    updategamma.rows(Extremevalue)=arma::ones(Extremevalue.n_rows)*std::pow(10,-10);
    
    arma::vec lastpar(beta.n_rows+1);
    lastpar.subvec(0,beta.n_rows-1)=beta;
    // lastpar(beta.n_rows)=std::log(theta);
    arma::vec updatePar(beta.n_rows+1);
    updatePar.subvec(0,beta.n_rows-1)=lastpar.subvec(0,beta.n_rows-1)-solve(SecondDeriv,FirstDeriv);
    updatePar(beta.n_rows)=theta;
    
    
    arma::field<arma::mat> result(2);
    result(0)=updategamma;
    result(1)=updatePar;
    return result;
    
}



// [[Rcpp::export]]
double LogProfileLikeli(const arma::uword&Indicator,const arma::mat&Data,const arma::mat&rules,const double&Tol,const arma::mat&beta,const double&theta,const arma::mat&gamma){
    if(Indicator<=beta.n_rows){
        int betadim=Data.n_cols-5;
        
        arma::mat tktemp=sort(unique(join_cols(Data.col(1),Data.col(2))));
        arma::vec tk=tktemp.col(0);
        tk=tk.subvec(1,tk.n_elem-2);
        
        double maxL=Data.col(1).max();
        arma::mat Rinspection=Data.col(2);
        Rinspection=Rinspection.rows(find(Rinspection<std::pow(10,10)));
        double maxR=Rinspection.max();
        arma::uvec tkfinity=find(tk<=(maxL)&&tk<=(maxR));
        tk=tk(tkfinity);
        
        arma::mat L=Data.col(1);
        arma::mat R=Data.col(2);
        
        arma::umat tL=TmatL(join_rows(Data.col(1),Data.col(2)),tk);
        arma::umat tR=TmatR(join_rows(Data.col(1),Data.col(2)),tk);
        arma::umat tLR=TmatLR(join_rows(Data.col(1),Data.col(2)),tk);
        arma::mat beta0=beta;
        arma::mat gamma0=gamma;
        double theta0=theta;
        double absdiff=1000;
        arma::mat betatheta;
        // arma::mat gammaold=gamma0;
        // arma::mat betaold=beta0;
        // double thetaold=theta0;
        // arma::mat Totalold=join_cols(gammaold,betaold,arma::ones(1,1)*thetaold);
        arma::field<arma::mat> newresult;
        int iter=1;
        do{
            arma::mat oldbeta=beta0;
            double oldtheta=theta0;
            arma::mat oldgamma=gamma0;
            newresult=UpdateonceProfileBeta(Indicator,rules,beta0,theta0,gamma0,Data,tL,tR,tLR);
            gamma0=newresult(0);
            betatheta=newresult(1);
            beta0=betatheta.submat(0,0,betadim-1,0);
            theta0=betatheta(betadim,0);
            arma::mat Totalnew=newresult(1);
            arma::mat Totalold=join_cols(oldbeta,arma::ones(1,1)*oldtheta);
            absdiff=accu(abs((Totalold-Totalnew)%pow(Totalold,-1)));
            iter=iter+1;
            // std::cout<<iter;
            // std::cout<<absdiff;
            
        } while (iter<500&absdiff>Tol);
        arma::field<arma::mat> finalresult(4);
        finalresult(0)=newresult(0);
        finalresult(1)=newresult(1);
        finalresult(2)=arma::ones(1,1)*absdiff;
        finalresult(3)=arma::ones(1,1)*iter;
        // return finalresult;
        double LogLike=LogLikeliValue(rules,beta0,theta0,finalresult(0),Data,tL,tR);
        return LogLike;
    }
    else{
        int betadim=Data.n_cols-5;
        
        arma::mat tktemp=sort(unique(join_cols(Data.col(1),Data.col(2))));
        arma::vec tk=tktemp.col(0);
        tk=tk.subvec(1,tk.n_elem-2);
        
        double maxL=Data.col(1).max();
        arma::mat Rinspection=Data.col(2);
        Rinspection=Rinspection.rows(find(Rinspection<std::pow(10,10)));
        double maxR=Rinspection.max();
        arma::uvec tkfinity=find(tk<=(maxL)&&tk<=(maxR));
        tk=tk(tkfinity);
        
        arma::mat L=Data.col(1);
        arma::mat R=Data.col(2);
        
        arma::umat tL=TmatL(join_rows(Data.col(1),Data.col(2)),tk);
        arma::umat tR=TmatR(join_rows(Data.col(1),Data.col(2)),tk);
        arma::umat tLR=TmatLR(join_rows(Data.col(1),Data.col(2)),tk);
        arma::mat beta0=beta;
        arma::mat gamma0=gamma;
        double theta0=theta;
        double absdiff=1000;
        arma::mat betatheta;
        // arma::mat gammaold=gamma0;
        // arma::mat betaold=beta0;
        // double thetaold=theta0;
        // arma::mat Totalold=join_cols(gammaold,betaold,arma::ones(1,1)*thetaold);
        arma::field<arma::mat> newresult;
        int iter=1;
        do{
            arma::mat oldbeta=beta0;
            double oldtheta=theta0;
            arma::mat oldgamma=gamma0;
            newresult=UpdateonceProfileTheta(rules,beta0,theta0,gamma0,Data,tL,tR,tLR);
            gamma0=newresult(0);
            betatheta=newresult(1);
            beta0=betatheta.submat(0,0,betadim-1,0);
            theta0=betatheta(betadim,0);
            arma::mat Totalnew=newresult(1);
            arma::mat Totalold=join_cols(oldbeta,arma::ones(1,1)*oldtheta);
            absdiff=accu(abs((Totalold-Totalnew)%pow(Totalold,-1)));
            iter=iter+1;
            // std::cout<<iter;
            // std::cout<<absdiff;
            
        } while (iter<500&absdiff>Tol);
        arma::field<arma::mat> finalresult(4);
        finalresult(0)=newresult(0);
        finalresult(1)=newresult(1);
        finalresult(2)=arma::ones(1,1)*absdiff;
        finalresult(3)=arma::ones(1,1)*iter;
        double LogLike=LogLikeliValue(rules,beta0,theta0,finalresult(0),Data,tL,tR);
        return LogLike;
        
    }
    
}


// [[Rcpp::export]]
double LogProfileLikeliInd(const arma::uword&Indicator,const arma::mat&Data,const arma::mat&rules,const double&Tol,const arma::mat&beta,const double&theta,const arma::mat&gamma){
    
        int betadim=Data.n_cols-5;
        
        arma::mat tktemp=sort(unique(join_cols(Data.col(1),Data.col(2))));
        arma::vec tk=tktemp.col(0);
        tk=tk.subvec(1,tk.n_elem-2);
        
        double maxL=Data.col(1).max();
        arma::mat Rinspection=Data.col(2);
        Rinspection=Rinspection.rows(find(Rinspection<std::pow(10,10)));
        double maxR=Rinspection.max();
        arma::uvec tkfinity=find(tk<=(maxL)&&tk<=(maxR));
        tk=tk(tkfinity);
        
        arma::mat L=Data.col(1);
        arma::mat R=Data.col(2);
        
        arma::umat tL=TmatL(join_rows(Data.col(1),Data.col(2)),tk);
        arma::umat tR=TmatR(join_rows(Data.col(1),Data.col(2)),tk);
        arma::umat tLR=TmatLR(join_rows(Data.col(1),Data.col(2)),tk);
        arma::mat beta0=beta;
        arma::mat gamma0=gamma;
        double theta0=theta;
        double absdiff=1000;
        arma::mat betatheta;
        // arma::mat gammaold=gamma0;
        // arma::mat betaold=beta0;
        // double thetaold=theta0;
        // arma::mat Totalold=join_cols(gammaold,betaold,arma::ones(1,1)*thetaold);
        arma::field<arma::mat> newresult;
        int iter=1;
        do{
            arma::mat oldbeta=beta0;
            arma::mat oldgamma=gamma0;
            
            newresult=UpdateonceProfileBetaInd(Indicator,rules,beta0,theta0,gamma0,Data,tL,tR,tLR);
            gamma0=newresult(0);
            betatheta=newresult(1);
            beta0=betatheta.submat(0,0,betadim-1,0);
            arma::mat Totalnew=newresult(1);
            arma::mat Totalold=oldbeta;
            absdiff=accu(abs((Totalold-Totalnew)%pow(Totalold,-1)));
            iter=iter+1;
            // std::cout<<iter;
            // std::cout<<absdiff;
            

            
        } while (iter<500&absdiff>Tol);
        arma::field<arma::mat> finalresult(4);
        finalresult(0)=newresult(0);
        finalresult(1)=newresult(1);
        finalresult(2)=arma::ones(1,1)*absdiff;
        finalresult(3)=arma::ones(1,1)*iter;
        // return finalresult;
        double LogLike=LogLikeliValue(rules,beta0,theta0,finalresult(0),Data,tL,tR);
        return LogLike;
    
    
}


// [[Rcpp::export]]
double FindSolution(const arma::uword&Indicator,const arma::mat&Data,const arma::mat&rules,const double&Tol,const arma::mat&beta,const double&theta,const arma::mat&gamma,
                    const double&Tol2,const double&TheConst,const arma::vec&Initial){
    if(Indicator<=beta.n_rows){
        double Initial_L=Initial(0);
        double Initial_R=Initial(1);
        double Initial_Ave;
        arma::mat beta0=beta;
        arma::mat gamma0=gamma;
        double theta0=theta;
        do{
            Initial_Ave=(Initial_L+Initial_R)/2;
            beta0(Indicator-1,0)=Initial_Ave;
            double plAve;
            plAve=LogProfileLikeli(Indicator,Data,rules,Tol,beta0,theta0,gamma0)-TheConst;
            if(plAve<0){
                Initial_L=Initial_Ave;
            }
            else{
                Initial_R=Initial_Ave;
            }
            
        } while (std::abs(Initial_R-Initial_L)>Tol2);
        return Initial_Ave;
    }
    else{
        double Initial_L=Initial(0);
        double Initial_R=Initial(1);
        double Initial_Ave;
        arma::mat beta0=beta;
        arma::mat gamma0=gamma;
        double theta0=theta;
        do{
            Initial_Ave=(Initial_L+Initial_R)/2;
            theta0=Initial_Ave;
            double plAve=LogProfileLikeli(Indicator,Data,rules,Tol,beta0,theta0,gamma0)-TheConst;
            if(plAve<0){
                Initial_L=Initial_Ave;
            }
            else{
                Initial_R=Initial_Ave;
            }
            
        } while (std::abs(Initial_R-Initial_L)>Tol2);
        return Initial_Ave;
    }
}


// [[Rcpp::export]]
double FindSolutionInd(const arma::uword&Indicator,const arma::mat&Data,const arma::mat&rules,const double&Tol,const arma::mat&beta,const double&theta,const arma::mat&gamma,
                    const double&Tol2,const double&TheConst,const arma::vec&Initial){
    
        double Initial_L=Initial(0);
        double Initial_R=Initial(1);
        double Initial_Ave;
        arma::mat beta0=beta;
        arma::mat gamma0=gamma;
        double theta0=theta;
        do{
            Initial_Ave=(Initial_L+Initial_R)/2;
            beta0(Indicator-1,0)=Initial_Ave;
            double plAve;
            
            plAve=LogProfileLikeliInd(Indicator,Data,rules,Tol,beta0,theta0,gamma0)-TheConst;
            if(plAve<0){
                Initial_L=Initial_Ave;
            }
            else{
                Initial_R=Initial_Ave;
            }
            
        } while (std::abs(Initial_R-Initial_L)>Tol2);
        return Initial_Ave;
    

}


// [[Rcpp::export]]
arma::mat ComputeCI(const arma::mat&Data,const arma::mat&rules,const double&Tol,const arma::mat&beta,const double&theta,const arma::mat&gamma,
                    const double&Tol2,const double&TheConst){
    arma::mat TotalEst=join_cols(beta,arma::ones(1,1)*theta);
    int Totaldim=beta.n_rows+1;
    arma::mat result=arma::zeros(Totaldim,2);
    arma::mat InL=arma::zeros(Totaldim,2);
    arma::mat InR=arma::zeros(Totaldim,2);
    InL.col(1)=TotalEst;
    InR.col(1)=TotalEst;
    InL.col(0)=InL.col(1)-2;
    InR.col(0)=InR.col(1)+2;
    InL(beta.n_rows,0)=0;
    for(arma::uword i=1;i<Totaldim+1;i++){
        result(i-1,0)=FindSolution(i,Data,rules,Tol,beta,theta,gamma,Tol2,TheConst,trans(InL.row(i-1)));
        result(i-1,1)=FindSolution(i,Data,rules,Tol,beta,theta,gamma,Tol2,TheConst,trans(InR.row(i-1)));
    }
    return result;
}



// [[Rcpp::export]]
arma::mat ComputeCIInd(const arma::mat&Data,const arma::mat&rules,const double&Tol,const arma::mat&beta,const double&theta,const arma::mat&gamma,
                    const double&Tol2,const double&TheConst){
    arma::mat TotalEst=(beta);
    int Totaldim=beta.n_rows;
    arma::mat result=arma::zeros(Totaldim,2);
    arma::mat InL=arma::zeros(Totaldim,2);
    arma::mat InR=arma::zeros(Totaldim,2);
    InL.col(1)=TotalEst;
    InR.col(1)=TotalEst;
    InL.col(0)=InL.col(1)-2;
    InR.col(0)=InR.col(1)+2;
    for(arma::uword i=1;i<Totaldim+1;i++){
        result(i-1,0)=FindSolution(i,Data,rules,Tol,beta,theta,gamma,Tol2,TheConst,trans(InL.row(i-1)));
        result(i-1,1)=FindSolution(i,Data,rules,Tol,beta,theta,gamma,Tol2,TheConst,trans(InR.row(i-1)));
    }
    return result;
}





