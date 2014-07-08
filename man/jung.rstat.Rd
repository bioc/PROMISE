\name{jung.rstat}
\alias{jung.rstat}
\Rdversion{3.1.0}
\title{Function to Compute Jung's Statistics }
\description{
  Compute statistic that measures the correlation of many continuous variables with a censored time-to-event variable
}
\usage{
jung.rstat(x, time.cens, strat = NULL)
}

\arguments{
  \item{x}{ 
  a data frame with row corresponding to probe set and column corresponding
  to subjects, the order of columns (subjects) should match the order of rows in \emph{time.cens}.
  }
  \item{time.cens}{ 
  a data frame with number of row equal to number of column in \emph{x}.
  It contains two columns with first for time and second for censor
  (1 = event, 0 = censored).   
  }
  \item{strat}{ 
  a vector of stratum to calculate stratified r-type association statistics, default = NULL.
  }
}

\value{
Returns a vector of Jun's r-type association statistics.
}
\references{ 
 Jung SH, Owzar K, and George SL (2005) A multiple testing procedure to 
 associate gene expression levels with survival. Stat Med 24:3077-88
 }
\author{ 
Stan Pounds \email{stanley.pounds@stjude.org}; Xueyuan Cao \email{xueyuan.cao@stjude.org}
}
\note{
The order of subjects in \emph{x} (column), \emph{time.cens}, and \emph{strat} should all match. 
The original statistic proposed by Jung, Owzar, and George can be written as a dot-product.  
The statistic returned by this routine is expressed in the form of a correlation statistic by dividing the dot product by the square root of the lengths of the two vectors in the numerator.
}
\seealso{
\code{\link{PROMISE}}
}
\examples{
## load sampExprSet.
data(sampExprSet)

## extract expression matrix from sampExprSet
Y <- exprs(sampExprSet)

## extract end point data from sampExprSet
time.cens <-pData(phenoData(sampExprSet))[, 3:4] 
strat <- pData(phenoData(sampExprSet))$strat

## compute Jung's r-type association statistics
jungstat <- jung.rstat(Y, time.cens, strat = strat)
}
\keyword{survival}

