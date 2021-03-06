% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_coefs.R
\name{test_coefs}
\alias{test_coefs}
\title{Standalone function to compare local and external datasets with a likelihood ratio test (LRT)}
\usage{
test_coefs(datasets, family = c("linear", "logistic", "cox"))
}
\arguments{
\item{datasets}{A list of datasets (each dataset having design matrix "x", response "y",
and potentially censoring vector "censor").
It is assumed the first indexed dataset is the local dataset.}

\item{family}{The model to perform an LRT for, one of "linear", "logistic", or "cox".}
}
\value{
A list of p-values corresponding to testing the null hypothesis that one shared set of parameters
is sufficient to describe the local dataset and each external dataset
}
\description{
Standalone function to compare local and external datasets with a likelihood ratio test (LRT)
}
