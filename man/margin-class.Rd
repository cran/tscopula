% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/margins.R
\docType{class}
\name{margin-class}
\alias{margin-class}
\alias{coef,margin-method}
\alias{sim,margin-method}
\alias{show,margin-method}
\title{Marginal model for time series}
\usage{
\S4method{coef}{margin}(object)

\S4method{sim}{margin}(object, n = 1000)

\S4method{show}{margin}(object)
}
\arguments{
\item{object}{an object of the class.}

\item{n}{length of realization.}
}
\description{
Class of objects for marginal models for stationary time series. The
object is given a name and there must exist functions pname, qname,
dname and rname. As well as the parameters of the distribution, dname must have the
logical argument log specifying whether log density should be computed.
}
\section{Methods (by generic)}{
\itemize{
\item \code{coef(margin)}: Coef method for margin class

\item \code{sim(margin)}: Simulation method for margin class

\item \code{show(margin)}: Show method for margin class

}}
\section{Slots}{

\describe{
\item{\code{name}}{name of the marginal model class.}

\item{\code{pars}}{a numeric vector containing the named parameters of the distribution
which are passed as arguments to pname, qname, dname and rname.}
}}

\examples{
new("margin", name = "gauss", pars = c(mu = 0, sigma = 1))
margmod <- margin("gauss", pars = c(mu = 0, sigma = 1))
sim(margmod, n = 500)
}
