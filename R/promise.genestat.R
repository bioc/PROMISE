`promise.genestat` <-
function(Y,ph.data,ph.pattern, strat=NULL, proj0=FALSE)

{
   m<-dim(Y)[1]
   res<-0
   nph<-dim(ph.pattern)[1]
   sum.wgt<-0
   ph.stats<-matrix(NA,m,nph)
   #ph.pvals<-matrix(NA,m,nph)
   stat.names<-pval.names<-rep("",nph)
   if (!proj0)
   {
     for (i in 1:nph)
     {
       ph.var<-unlist(strsplit(as.character(ph.pattern$endpt.vars[i]),split=","))
       temp.ph<-ph.data[,ph.var]

       ph.stat<-as.character(ph.pattern$stat.func[i])
       ph.wgt<-ph.pattern$stat.coef[i]
       ph.call<-call(ph.stat,Y,temp.ph, strat)
       ph.val<-eval(ph.call)
       ph.stats[,i]<-ph.val
       stat.names[i]<-paste(unlist(c(unlist(ph.var),"stat")),collapse=".")
       res<-res+ph.wgt*ph.val
       sum.wgt<-sum.wgt+abs(ph.wgt)
     }
     res<-res/sum.wgt
   }
   if (proj0)
   { 
     for (i in 1:nph)
     {
       ph.var<-unlist(strsplit(as.character(ph.pattern$endpt.vars[i]),split=","))
       temp.ph<-ph.data[,ph.var]

       ph.stat<-as.character(ph.pattern$stat.func[i])
       ph.call<-call(ph.stat,Y,temp.ph, strat)
       ph.val<-eval(ph.call)
       ph.stats[,i]<-ph.val
       stat.names[i]<-paste(unlist(c(unlist(ph.var),"stat")),collapse=".")
       res<-res+ph.val^2
     }    
   }  
   colnames(ph.stats)<-stat.names
   stat<-cbind(ph.stats,PROMISE.stat=res)
   return(stat)
}

