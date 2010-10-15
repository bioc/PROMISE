`avg.abs.genestat` <-
function(gene.res,probes,GS.data)
#gene.res: individual gene result
#probes: probeid
#GS.data: gene set definition data with probe and set variable
{
    ugs<-unique(GS.data[,2])
    ngs<-length(ugs)
    nph<-dim(gene.res)[2]
    res<-matrix(NA,ngs,nph)
    A<-abs(gene.res)
    for (i in 1:ngs)
    {
       gs<-is.element(GS.data[,2],ugs[i])
       in.gs<-is.element(probes,GS.data[gs,1])
       if (sum(in.gs)==1) res[i,] <-mean(A[in.gs,])
       else  res[i,]<-colMeans(A[in.gs,])
    }
    colnames(res)<-colnames(A)
    return(res)
}

