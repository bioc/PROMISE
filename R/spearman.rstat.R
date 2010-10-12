`spearman.rstat` <-
function(Y,  x,  strat=NULL)
# Y: Array data matrix, subjects in columns, variables in rows
# x: phenotype vector, entries for subjects, same order as columns of Y
# strat=NULL: stratification vector, entries for subjects, same order as columns of Y
{
  if (is.null(strat)) strat<-rep(1, length(x))
  # Some initialization
	ustrat<-unique(strat)
	nstrat<-length(ustrat)
	stat.sum<-0
	sum.n<-0

  # Loop over stratum
	for (i in 1:nstrat)
	{
    this.strat<-(strat==ustrat[i])
		this.x<-x[this.strat]
		this.Y<-Y[,this.strat]
    this.n<-sum(this.strat)

       this.x<-unlist(this.x)
       this.x.miss<-is.na(this.x)
       this.x<-this.x[!this.x.miss]
       this.Y<-this.Y[,!this.x.miss]

       r<-rank(this.x)
       r.mn<-mean(r)
       r.sd<-sd(r)
       r.cnt<-(r-r.mn)/r.sd
       r.cnt[r.sd==0]<-0

       # Now identify missing Y
       this.Y.miss<-is.na(this.Y)

       # Now rank Y
       tR<-apply(this.Y,1,rank)
       R<-t(tR)
       R[this.Y.miss]<-NA
       n<-rowSums(!this.Y.miss)
       R.mn<-rowMeans(R,na.rm=TRUE)
       R.dev<-(R-R.mn)
       R.sd<-sqrt(rowMeans(R.dev^2,na.rm=TRUE))
       R.cnt<-(R-R.mn)/R.sd
       R.cnt[this.Y.miss]<-0
       R.cnt[R.sd==0]<-0
       res<-(R.cnt%*%r.cnt)/n
       res<-as.vector(res)
       
		stat.sum<-stat.sum+this.n*res
		sum.n<-sum.n+this.n
	}
  stat.sum<-as.vector(stat.sum/sum.n)
	return(stat.sum)
}

