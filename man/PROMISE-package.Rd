\name{PROMISE-package}
\alias{PROMISE-package}
\Rdversion{2.9.0}
\docType{PROMISE}
\title{
PRojection Onto the Most Interesting Statistical Evidence
}
\description{
a tool to identify genomic geatures with a
specific biologically interesting pattern of associations
with multiple endpoint variables
}
\details{
\tabular{ll}{
Package: \tab PROMISE\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2010-1-12\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
The PROMISE (PRojection Onto the Most Interesting Statistical Evidence) is performed by calling function PROMISE.
The array data and endpoint data are passed through an ExpressionSet; the gene set definition is passed throught a GeneSetCollection, 
and PROMISE definition is passed through a data frame. \emph{promise.genestat} and \emph{avg.abs.genestat} are called internally by PROMISE.
Two R routines for calculating association statistics with individual endpoint variable(\emph{jung.rstat} and \emph{spearman.rstat}) are 
provided in this version. Users could provide their own R routines written in a similar fashion. 
}
\author{
Stan Pounds \email{stanley.pounds@stjude.org}; Xueyuan Cao \email{xueyuan.cao@stjude.org}

Maintainer: Stan Pound \email{stanley.pounds@stjude.org}; Xueyuan Cao \email{xueyuan.cao@stjude.org}
}
\references{
Jung, S-H, Owzar K, and Goerge SL (2005) A multiple testing procedure to associate gene expression levels with survival.  
Biostatistics 24: 3077-3088.

Goeman JJ  and Buhlmann P (2007) Analyzing gene expression data in terms of gene sets: methodological issues. 
Bioinformatics  23: 980-987.

Pounds S, Cheng C, Cao X, Crews KR, Plunkett W, Gandhi V, Rubnitz J, Ribeiro RC, Downing JR, and Lamba J (2009) 
PROMISE: a tool to identify genomic features with a specific biologically interesting pattern of associations 
with multiple endpoint variables. Bioinformatics 25: 2013-2019
}
\keyword{ package }
\examples{
## load sampExprSet, sampGeneSet, phPatt.
data(sampExprSet)
data(sampGeneSet)
data(phPatt)

## Perform PROMISE procedure without GSEA
test1<-PROMISE(exprSet=sampExprSet,
             geneSet=NULL, 
             promise.pattern=phPatt, 
             strat.var=NULL, 
             seed=13, 
             nperms=10)
             
## Perform PROMISE procedure with GSEA             
res<-PROMISE(exprSet=sampExprSet,
             geneSet=sampGeneSet, 
             promise.pattern=phPatt, 
             strat.var=NULL, 
             seed=13, 
             nperms=10)
}

