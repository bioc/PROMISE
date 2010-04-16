\name{spearman.rstat}
\Rdversion{2.9.0}
\alias{spearman.rstat}
\title{
Function to Calculate Spearman Correlation Statistics
}
\description{
A function to calculate Spearman rank correlation of each gene in an array data with a continuous variable
}
\usage{
spearman.rstat(Y, x, strat = NULL)
}
\arguments{
  \item{Y}{
a numeric data frame. Each row gives values of one genomic variable.
}
  \item{x}{
a vector of continuous variable.
}
  \item{strat}{
a vector of stratum to calculate stratified correlation statistics, default = NULL.
}
}
\value{
Return a vector of spearman rank correlation statistics.
}
\references{
Spearman C. (1904) The proof and measurement of association between two things. Amer. J. Psychol. 15: 72-101
}
\author{   
Stan Pounds \email{stanley.pounds@stjude.org}; Xueyuan Cao \email{xueyuan.cao@stjude.org}
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
x <- pData(phenoData(sampExprSet))$drugLevel
strat <- pData(phenoData(sampExprSet))$strat

## Calculte spearman correlation statistics
test <- spearman.rstat(Y, x, strat = strat)
}
\keyword{univar}
