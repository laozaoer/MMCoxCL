\name{SIM_DATA}
\alias{SIM_DATA}
\title{
Generating datasets for the main function Cox_est
}
\description{
SIM_DATA is used to simulate clustered or non-clustered interval-censored data. This function makes use of the Cox model to generate the time-to-event. 
}
\usage{
SIM_DATA(beta,  theta, n, H, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{beta}{
A \eqn{p\times 1} vector representing the regression parameter for  covariates.
}

  \item{theta}{
A positive scalar representing the standard error of the clustering effect. When it is zero, the subjects are all independent (no clustering). 
}
  \item{n}{
An integer that denotes the number of units (clusters).
}
  \item{H}{
A function representing the \eqn{H} function included the generalized odds rate model. \eqn{H} should be a non-negative and non-decreasing function with \eqn{H(0)=0}.
}

 \item{k}{
In this function, the size of each unit is generated from a discrete uniform distribution from 1 to \eqn{k}. Users should input the specific \eqn{k}.

 }
}

\value{
Function SIM_DATA returns a matrix, where the first column represents the unit index of each observation; the next two columns represent the left and right inspection times, respectively; the next two column represent the left and right censoring indicator, respectively; the rest of the columns represent the covariates. The first covariate is generated from Uniform\eqn{(-1,1)} and the second covariate is generated from Bernoulli(0.5). 



}