\name{promise.genestat}
\alias{promise.genestat}
\Rdversion{2.9.0}
\title{
Function to Calculate PROMISE Statistics
}
\description{
a function to calculate individual gene and PROMISE statistics for a defined pattern of association
}
\usage{
promise.genestat(Y, ph.data, ph.pattern, strat = NULL)
}
\arguments{
  \item{Y}{
  a data frame with row corresponding to probe set and column corresponding
  to subjects, the order of column should match order of row in \emph{ph.data}.
  }
  \item{ph.data}{
  a data frame with rows corresponding to subjects and columns 
  corresponding to endpoint variables.
  }
  \item{ph.pattern}{
  a data frame with column headers: \emph{stat.coef}, \emph{stat.func}, \emph{endpt.vars}. The \emph{stat.coeff} column gives 
  the coefficients for combining the statistics of association of genomic variable with individual endpoint 
  variable into the ultimate PROMISE statistic. The \emph{stat.func} column gives the name of the R routine that 
  computes the test statistic of association of the end point variables.
  \emph{jung.rstat} and \emph{spearman.rstat} are provided. Users can provide their own routines accordingly. 
  The \emph{endpt.vars} column gives the name(s) of variable(s) in  
  \emph{ph.data} needed to compute each term of the PROMISE statistic. A \bold{comma} without a space should be used 
  to separate multiple variables that correspond to the same term in the association pattern definition.       
  }
  \item{strat}{
  a vector of stratum to calculate stratified statistics. The default is NULL.
  }
}

\value{
 a matrix of statistics. Each row gives gene's statistics of each individual endpoint and the PROMISE analysis defined in \emph{ph.pattern}.
}
\references{
Pounds S, Cheng C, Cao X, Crews KR, Plunkett W, Gandhi V, Rubnitz J, Ribeiro RC, Downing JR, and Lamba J (2009) 
PROMISE: a tool to identify genomic features with a specific biologically interesting pattern of associations 
with multiple endpoint variables. Bioinformatics 25: 2013-2019
}
\author{
Stan Pounds \email{stanley.pounds@stjude.org}; Xueyuan Cao \email{xueyuan.cao@stjude.org}
}
\note{
a function internally called by PROMISE.
}
\seealso{\code{\link{PROMISE}}}
\examples{
## load sampExprSet, phPatt.
data(sampExprSet)
data(phPatt)

Y <- exprs(sampExprSet)
ph.data <- pData(phenoData(sampExprSet))

test <- promise.genestat(Y, ph.data, phPatt, strat=ph.data[, 5])
}
\keyword{multivariate}

