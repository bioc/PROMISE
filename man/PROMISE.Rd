\name{PROMISE}
\alias{PROMISE}
\Rdversion{3.1.0}
\title{
PRojection onto the Most Interesting Statistical Evidence
}
\description{
Perform permutation-based test to identify genes with expression levels having a specific biologically interesting pattern of associations
with multiple endpoint variables
}
\usage{
PROMISE(exprSet, geneSet=NULL, promise.pattern, strat.var=NULL, proj0=FALSE, seed=13, nbperm=FALSE, max.ntail=100,
        nperms=10000)
}
\arguments{
  \item{exprSet}{
   an ExpressionSet class contains minimum of \emph{exprs} (\emph{expression matrix}) and \emph{phenoData} (\emph{AnnotatedDataFrame of end point data}). 
   Please refer to \emph{Biobase} for details on how to create such an ExpressionSet.}
  \item{geneSet}{
   a GeneSetCollection class with minimum of setName and geneIDs for each GeneSet. 
   Please refer to \emph{GSEABase} for how to create such a GeneSetCollection class. The default is NULL which will perform no gene set enrichment analysis.}  
  \item{promise.pattern}{
  a data frame defining the association pattern of interest. The column names must be \emph{stat.coef}, \emph{stat.func}, and \emph{endpt.vars}.
  The \emph{stat.coef} column gives the coefficients for combining the statistics of association of genomic variables with individual endpoint variable into the ultimate PROMISE statistic. If \emph{proj0=TRUE}, the \emph{stat.coeff} is ignored.  
  The \emph{stat.func} column gives the name of the R routine that computes the test statistic of association of the endpoint variables.
  Two R routines (\emph{jung.rstat} and \emph{spearman.rstat})are provided. Users can provide their own routine accordingly. 
  The \emph{endpt.vars} column gives the name(s) of variable(s) in the endpoint data file needed to compute each term of 
  the PROMISE statistic. A common without a space should be used to separate multiple variables that correspond to the same term in the association pattern definition.
  } 
  \item{strat.var}{
  the name or numeric value of stratum variable in \emph{exprSet} for stratified analysis. The default is NULL which performs an unstratified analysis.
  }
  \item{proj0}{
  indicator of whether projection to 0 is performed. It takes two valid values: TRUE or FALSE. If \emph{proj0=TRUE}, PROMISE statistics is the sum of squares of individual statitics and the \emph{stat.coeff} in \emph{promise.pattern} is ignored. The default is FALSE.
  }
  \item{seed}{
  initial seed of random number generator. The default is 13.}
  \item{nbperm}{
  indicator of fast permuation using negative binomial strategy, taking two valid values: FALSE or TRUE. The default is FALSE.}
  \item{max.ntail}{
  number of sucess if \emph{nbperm = T}. Further permutation will not be performed for gene(s) or gene set(s) which \emph{max.ntail} permutated statistics are greater or equal to the observed statistics, The default is 100}
  \item{nperms}{
  number of permutations. The default is 10,000.}
}
\value{
  \item{$generes }{individual genes' test statistics and p-values for each individual endpoint 
  and PROMISE analysis. If \emph{nbperm=T}, the last column contains number of permuations for each gene.
  }
  \item{$setres}{gene set's test statistics and p-values for each individual endpoint 
  and PROMISE analysis. If \emph{nbperm=T}, the last column contains number of permuations for gene set. If \emph{geneSet} is NULL, the value of this component is also \emph{NULL}. 
  }
}
\references{
Pounds S, Cheng C, Cao X, Crews KR, Plunkett W, Gandhi V, Rubnitz J, Ribeiro RC, Downing JR, and Lamba J (2009) 
PROMISE: a tool to identify genomic features with a specific biologically interesting pattern of associations 
with multiple endpoint variables. Bioinformatics 25: 2013-2019
}
\author{
Stan Pounds \email{stanley.pounds@stjude.org}; Xueyuan Cao \email{xueyuan.cao@stjude.org}
 }
\seealso{\code{\link{jung.rstat}} \code{\link{avg.abs.genestat}} 
\code{\link{promise.genestat}}  \code{\link{spearman.rstat}}
\code{\link{promise.pattern}}
}
\examples{
## load sampExprSet, sampGeneSet, phPatt.
data(sampExprSet)
data(sampGeneSet)
data(phPatt)

## Perform PROMISE procedure without GSEA
test1 <- PROMISE(exprSet=sampExprSet,
                 geneSet=NULL, 
                 promise.pattern=phPatt, 
                 strat.var=NULL,
                 proj0=FALSE,
                 nbperm=FALSE, 
                 max.ntail=10,
                 seed=13, 
                 nperms=100)
             
## Perform PROMISE procedure with GSEA and using fast permuation              
test2 <- PROMISE(exprSet=sampExprSet,
                 geneSet=sampGeneSet, 
                 promise.pattern=phPatt, 
                 strat.var=NULL,
                 proj0=TRUE,
                 nbperm=TRUE, 
                 max.ntail=10,
                 seed=13, 
                 nperms=100)
}        
\keyword{multivariate}

