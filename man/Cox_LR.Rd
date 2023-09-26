\name{Cox_LR}
\alias{Cox_LR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performing the likelihood ratio test for the presence of clustering effect.
}
\description{
Cox_LR is used to calculate the test statistic and p-value of the the likelihood ratio test for the presence of clustering effect.
}
\usage{
Cox_LR(Data)
}
\arguments{
  \item{Data}{
A data frame, where the first column represents the clustering  
index of each observation (level 2 index); the next two columns
represent the left and right inspection times, respectively; the next two
columns represent the left and right censoring indicator, respectively.
the rest of the columns represent the covariates. Just understand that if one of the covariates is a cluster-level covariate, then its value is the same for all level 1 units within any 
level 2 unit (cluster). If there are any NA values in the data frame, then the entire row where the NA occurs will be ignored in the analysis.    
}

}

\value{
Function Cox_est returns a dataframe, where the first column contains the test statistic and the second column contains the corresponding p-value.
}

\examples{
library(MMCoxCL)
set.seed(1)
H=function(t) t^3 # specify the H function



#Scenario 1
#There is clustering effect.
H=function(t) t^(3)
data=SIM_DATA(beta=c(0.5,0.5),theta=0.5,n=100,H=H,k=10)
#performe the likelihood ratio test
result=Cox_LR(data)


}