\name{Only_Cox_est}
\alias{Only_Cox_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating parameters of the Cox model with clustered interval-censored data.
}
\description{
Only_Cox_est is used to analyze  clustered interval-censored data by fitting the Cox model. This is the main function of this package. It is assumed that a clustered dataset has \eqn{n} units and there are \eqn{n_i} subjects within the \eqn{i}th unit.  
}
\usage{
Cox_est(Data, Cluster.ind = TRUE, Tol = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
A data frame, where the first column represents the clustering  
index of each observation (level 2 index); the next two columns
represent the left and right inspection times, respectively; the next two
columns represent the left and right censoring indicator, respectively.
the rest of the columns represent the covariates. Just understand that if one of the covariates is a cluster-level covariate, then its value is the same for all level 1 units within any 
level 2 unit (cluster). If there are any NA values in the data frame, then the entire row where the NA occurs will be ignored in the analysis.    
}
  \item{Cluster.ind}{
Logical. If TRUE, the parameters will be estimated along with clustering effect. If it is FALSE, then all subjects are treated as independent and there is no clustering effect.  The default value is TRUE.
}

  \item{Tol}{
This denotes the sum of the absolute relative differences of the estimates at two
successive iterations. It is used to check the convergence of the parameter estimates.
The default value is 1e-5.
}
}

\value{
Function Cox_est returns a dataframe, where the first column contains the estimate of the regression parameters.
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
H=function(t) t^(3)
data=SIM_DATA(beta=c(0.5,0.5),theta=0.5,n=100,H=H,k=10)
#estimate the parameters
result=Cox_est(data)

#Scenario 2
#There is not clustering effect.
H=function(t) t^(3)
data=SIM_DATA(beta=c(0.5,0.5),theta=0,n=100,H=H,k=10)
#estimate the parameters
result=Cox_est(data,Cluster.ind=FALSE)

}
