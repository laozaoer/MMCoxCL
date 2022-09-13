\name{Cox_est}
\alias{Cox_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating parameters of the Cox model with clustered interval-censored data.
}
\description{
Cox_est is used to analyze  clustered interval-censored data by fitting the Cox model. This is the main function of this package. It is assumed that a clustered dataset has \eqn{n} units and there are \eqn{n_i} subjects within the \eqn{i}th unit.  
}
\usage{
Cox_est(Data, Cluster.ind = TRUE, CIlevel = 0.95, Tol = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
A matrix, where the first column represents the unit index of each observation; the next two columns represent the left and right inspection times, respectively; the next two column represent the left and right censoring indicator, respectively; the rest of the columns represent the covariates. If one of the covariates is cluster-level covariate, just let this covariate value as same for all the observations within the same cluster.
}
  \item{Cluster.ind}{
Logical. If TRUE, the parameters will be estimated along with clustering effect. If it is FALSE, then all subjects are treated as independent and there is no clustering effect.  The default value is TRUE.
}
  \item{CIlevel}{
A scalar containing the confidence level of the confidence interval. The defaul value is 0.95.
}
  \item{Tol}{
This denotes the sum of the absolute relative differences of the estimates at two
successive iterations. It is used to check the convergence of the parameter estimates.
The default value is 1e-5.
}
}

\value{
Function Cox_est returns a dataframe, where the first column contains the estimate of the regression parameters and the next two columns contain the lower and upper bounds of the confidence interval, respectively.
}
\references{
%% ~put references to the literature/web site here ~
}

\examples{
library(MMCoxCL)
set.seed(1)
H=function(t) t^3 # specify the H function



#Scenario 1
#There is clustering effect.
data=SIM_DATA(beta=c(0.5,0.5),theta=0.5,n=100,H=H,k=10)
#estimate the parameters
result=Cox_est(data)

#Scenario 2
#There is not clustering effect.
data=SIM_DATA(beta=c(0.5,0.5),theta=0,n=100,H=H,k=10)
#estimate the parameters
result=Cox_est(data,Cluster.ind=FALSE)

}