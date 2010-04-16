\name{avg.abs.genestat}
\alias{avg.abs.genestat}
\Rdversion{2.9.0}
\title{
Function to Compute Gene Set Statistics
}
\description{
A function to calculate the mean of absolute values of statistics based on a gene set definition
}
\usage{
avg.abs.genestat(gene.res, probes, GS.data)
}
\arguments{
  \item{gene.res}{
  a data frame. Each row gives test statistics for a genomic variable. Each column corresponds to
  an endpoint variable.
  }
  \item{probes}{
  a vector that links the \emph{gene.res} to \emph{GS.data}.                                         
  }
  \item{GS.data}{
  a data frame with first column for probe set identifier and second column for gene set identifier. 
  Each row assigns a probe set to a gene set. Each probe set may be assigned to multiple gene sets or no gene set at all.
  }
}
\value{
 Return a matrix of statistics. Each row gives the mean absolute value of test statistics of genes belonging to a gene set. 
 The columns are same as in \emph{gene.res}.
 }
\references{
Goeman JJ  and Buhlmann P (2007) Analyzing gene expression data in terms of gene sets: methodological issues. 
Bioinformatics  23: 980-987.
}
\author{
Stan Pounds \email{stanley.pounds@stjude.org}; Xueyuan Cao \email{xueyuan.cao@stjude.org}
 }
\seealso{
\code{\link{PROMISE}}
}
\note{
A function internally called by \emph{PROMISE}. 
}
\examples{

## load sampExprSet sampGeneSet.
data(sampExprSet)
data(sampGeneSet)

## extract expression matrix from sampExprSet
Y <- exprs(sampExprSet)
probes <- rownames(Y)

## convert sampGeneSet to a data frame
GS.data <- NULL
for (i in 1:length(sampGeneSet)){
    tt <- sampGeneSet[i][[1]]
    this.name <- unlist(geneIds(tt))
    this.set <- setName(tt)
    GS.data <- rbind.data.frame(GS.data,
                cbind.data.frame(featureID=as.character(this.name),
                                 setID=rep(as.character(this.set),
                                 length(this.name))))
}

## Calculate the mean of absolute values of statistics
## This is only a demo, probe expression values are used 
##in stead of statistics 
test <- avg.abs.genestat(Y, probes, GS.data)
}
\keyword{multivariate}
