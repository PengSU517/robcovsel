% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roblasso.R
\name{copulalasso}
\alias{copulalasso}
\title{Copula lasso}
\usage{
copulalasso(
  x,
  y,
  cor.method = "gaussrank",
  scale.method = "qn",
  center.method = "median",
  adaptive = TRUE
)
}
\arguments{
\item{x}{input design matrix}

\item{y}{input response vector}

\item{cor.method}{could be "pearson" or "gaussrank"}

\item{scale.method}{"sd" or "qn"}

\item{center.method}{mean or median}

\item{adaptive}{lasso or adaptive lasso}
}
\value{
betahat_opt, the optimal estimation of beta obtained from this algorithm

lambda_opt is the optimal tuning parameter and sigma_opt is the optimal estimation of sigma.

The output also includes the estimated correlation matrix, the estimated covariance matrix others.
}
\description{
Copula lasso
}
