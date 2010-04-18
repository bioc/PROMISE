\name{promise.pattern}
\alias{promise.pattern}
\title{
PROMISE pattern
}
\description{
PROMISE pattern is a data frame of association pattern definition, consisting of three columns. 
}
\format{
PROMISE pattern:  The column names must be \emph{stat.coef}, \emph{stat.func}, and \emph{endpt.vars}.
 
\emph{stat.coef} column gives the coefficients for combining the statistics of association of genomic variable with individual endpoint variable into the ultimate PROMISE statistic.
 
\emph{stat.func} column gives the name of the R routine that computes the test statistic of association of the endpoint variables.
  Two R routines (\emph{jung.rstat} and \emph{spearman.rstat})are provided in current release. Users can provide their own routine accordingly.
   
\emph{endpt.vars} column gives the name(s) of variable(s) in the endpoint data frame needed to compute each term of 
  the PROMISE statistic. If more than one variables involve in one term, they should be separated by a \bold{comma} without space. 
}

\references{
Pounds S, Cheng C, Cao X, Crews KR, Plunkett W, Gandhi V, Rubnitz J, Ribeiro RC, Downing JR, and Lamba J (2009) 
PROMISE: a tool to identify genomic features with a specific biologically interesting pattern of associations 
with multiple endpoint variables. Bioinformatics 25: 2013-2019
}
\author{
Stan Pounds \email{stanley.pounds@stjude.org}; Xueyuan Cao \email{xueyuan.cao@stjude.org}
 }
 \seealso{\code{\link{PROMISE}} 
}

\keyword{misc}
    
