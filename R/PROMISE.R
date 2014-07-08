`PROMISE` <-
function(exprSet, geneSet=NULL, promise.pattern, strat.var=NULL, proj0=FALSE, seed=13, nbperm=FALSE, max.ntail=100, nperms=10000)

{
     # extract data components from objects
     ph.pattern<-promise.pattern
     endpt.data<-pData(phenoData(exprSet))
     array.data<-exprs(exprSet)
     if (is.data.frame(array.data)) array.data<-as.matrix(array.data)
     if (!is.numeric(array.data))
     {
        print ("ERROR: The array data is not numeric")
       return()
     } 
     if (!is.null(geneSet)) {
          geneset.data<-NULL
          for (i in 1:length(geneSet)){
              tt<-geneSet[i][[1]]
              this.name<-unlist(geneIds(tt))
              this.set<-setName(tt)
              geneset.data<-rbind.data.frame(geneset.data, cbind.data.frame(featureID=as.character(this.name),
                                          setID=rep(as.character(this.set), length(this.name))))
          }
      }
     #geneset.data<-geneset.data[order(geneset.data[, 1]), ]
     # Match array data and endpttype data by subject identifiers
     array.names<-dimnames(array.data)[[2]]
     endpt.names<-dimnames(endpt.data)[[1]]

     array.inc<-is.element(array.names,endpt.names)
     endpt.inc<-is.element(endpt.names,array.names)
     nsamp.inc<-sum(array.inc)
     message(paste("No. of Matching array ids in endpt data and array data:",nsamp.inc))
     
     if (nsamp.inc<=4)
     {
        message(paste("Array IDs in endpt data:"))
        message(endpt.names)
        message("Array IDs in array data:")
        message(array.names)
        stop("Error: Unable to match array identifiers in endpt data and array data.")
     }
     
          
     # Subset on array data and phenotype data with matching subject IDs
     Y<-array.data[,array.inc]
     #Y<-Y[order(dimnames(Y)[[1]]), ]
     X<-endpt.data[endpt.inc,]
     array.names<-array.names[array.inc]
     endpt.names<-endpt.names[endpt.inc]

     array.ord<-order(array.names)
     endpt.ord<-order(endpt.names)

     Y<-Y[,array.ord]
     X<-X[endpt.ord,]
     Y<-as.matrix(Y)

     array.names<-array.names[array.ord]
     endpt.names<-endpt.names[endpt.ord]

     if (any(array.names!=endpt.names))
     {
        name.mtx<-cbind(endpt.arrayid=endpt.names,array.arrayid=array.names)
        message("----------------- IDs ------------------")
        message(name.mtx)
        stop("Error: Unable to align IDs in array data and endpt data")
     }
     
     endpt.endptcol<-is.element(dimnames(X)[[2]],unique(unlist(strsplit(c(as.character(ph.pattern$endpt.vars),
                                         strat.var), ","))) )
     phtype<-X[,endpt.endptcol]
     names(phtype)<-dimnames(endpt.data)[[2]][endpt.endptcol]
     if (!is.numeric(phtype[, dimnames(phtype)[[2]]!=strat.var])&!is.data.frame(phtype[, dimnames(phtype)[[2]]!=strat.var]))
     {
        print ("ERROR: The clinical data is not numeric or data frame")
       return()
     } 

     if (is.vector(phtype)) phtype<-matrix(phtype,length(phtype),1)
     
     # Now extract probe set names
     probes<-dimnames(array.data)[[1]]
     nprobes<-length(probes)
     
     if (nbperm)
     {      
       if (!is.null(geneSet)) 
       {

         #Now extract probe and gene set names
         gs.probes<-geneset.data[,1]
         gs.sets<-geneset.data[,2]

         uset<-unique(gs.sets)
         nset<-length(uset)
        
         # Now compute results for observed data
         message(paste("Computing Stats for Observed Data:",date()))
         message(paste("Computing gene-specific statistics:",date()))
         if (is.null(strat.var)) strat<-NULL  
         else  strat<-phtype[,strat.var] 
         gene.res<-promise.genestat(Y,phtype,ph.pattern, strat, proj0=proj0)

         m.nph<-dim(gene.res)
         gene.pvals<-matrix(0,m.nph[1],m.nph[2])
         
         message(paste("Computing gene-set enrichment statistics:",date()))
         enrich.res<-avg.abs.genestat(gene.res,probes, geneset.data)
     
         ngs.nph<-dim(enrich.res)
         enrich.pvals<-matrix(0,ngs.nph[1],ngs.nph[2])
     
         #set seed for sampling subjects
         message(paste("Set seed: seed = ",seed))
         set.seed(seed)
         
         probe.keep<-rep(T,nprobes)
         probe.done<-rep(F,nprobes)
         probe.ntail<-matrix(0,m.nph[1],m.nph[2])
         probe.apt.ntail<-rep(0,m.nph[1])
         
         set.keep<-rep(T,ngs.nph[1])
         set.done<-rep(F,ngs.nph[1])
         set.ntail<-matrix(0,ngs.nph[1],ngs.nph[2])
         set.apt.ntail<-rep(0,ngs.nph[1])

         set.probe.ntail<-rep(NA,nprobes)
         set.probe.keep<-rep(NA,nprobes)
         set.probe.keep[is.element(probes, geneset.data[, 1])]<-T
         set.probe.done<-rep(F,nprobes)
         
         num.perms<-rep(NA,nprobes)
         set.perms<-rep(NA,nset)
         
         i<-0
         while((i<nperms)&&(any(probe.keep)))
         {
           i<-i+1
           if (is.null(strat.var)) perm.ind<-sample(nsamp.inc,replace=FALSE)
           else
           {
             strat<-X[,strat.var]
             ustrat<-sort(unique(strat))
             nstrat<-length(ustrat)
             ind<-length(strat)
             for (j in 1:nstrat)
             {
               this.strat<-(strat == ustrat[j])
               perm.ind[this.strat]<-sample((1:ind)[this.strat],replace=FALSE)
             }
           }
           perm.phtype<-phtype[perm.ind,]
           
           if (is.null(strat.var))  perm.strat<-NULL 
           else  perm.strat<- perm.phtype[,strat.var]
                    
           gene.temp<-promise.genestat(Y[probe.keep, ],perm.phtype, ph.pattern, perm.strat, proj0=proj0)
           set.probe.keep2<-set.probe.keep[!is.na(set.probe.keep)]
           enrich.temp<-avg.abs.genestat(gene.temp,probes[probe.keep],geneset.data[set.probe.keep2, ])
           
           set.ntail[set.keep,]<-set.ntail[set.keep,] + (abs(enrich.temp)>=abs(enrich.res[set.keep, ]))
           set.apt.ntail[set.keep]<-set.ntail[set.keep, ngs.nph[2]]
           
           probe.ntail[probe.keep,]<-probe.ntail[probe.keep,]+(abs(gene.temp)>=abs(gene.res[probe.keep,]))
           probe.apt.ntail[probe.keep]<- probe.ntail[probe.keep, m.nph[2]]
           
           for (j in 1:sum(set.keep)){ 
             set.probe.ntail[!is.na(set.probe.keep)][is.element(probes[!is.na(set.probe.keep)], geneset.data[set.probe.keep2,][geneset.data[set.probe.keep2, 2]==uset[set.keep][j], 1]) ] <-
               set.apt.ntail[set.keep][j]  
           }
           
           probe.done[probe.keep]<-(probe.apt.ntail[probe.keep]>=max.ntail & set.probe.ntail[probe.keep] >= max.ntail) |  
                                   (probe.apt.ntail[probe.keep]>=max.ntail& is.na(set.probe.keep[probe.keep]))  
           set.probe.done[probe.keep & !is.na(set.probe.keep)]<-probe.done[probe.keep & !is.na(set.probe.keep)]
           
           for (j in 1:sum(set.keep)){
               set.done[set.keep][j] <- set.apt.ntail[set.keep][j]>=max.ntail & 
                 min(probe.apt.ntail[!is.na(set.probe.keep)][is.element(probes[!is.na(set.probe.keep)], geneset.data[set.probe.keep2,][geneset.data[set.probe.keep2, 2]== uset[set.keep][j], 1])])>=max.ntail
           }
           
           num.perms[probe.done&probe.keep]<-i
           set.perms[set.done&set.keep]<-i
           gene.pvals[probe.done&probe.keep, ]<-probe.ntail[probe.done&probe.keep,]/i

           enrich.pvals[set.done&set.keep, ] <-set.ntail[set.done&set.keep, ]/i  
           probe.keep[probe.done]<-F
           set.probe.keep[!is.na(set.probe.keep)& set.probe.done]<-F
           set.keep[set.done]<-F
         }
         
         if (any(probe.keep))
         {
           gene.pvals[probe.keep,]<-probe.ntail[probe.keep,]/nperms
           num.perms[probe.keep]<-nperms
         }
                  
         if (any(set.keep))
         {
           enrich.pvals[set.keep,]<-set.ntail[set.keep,]/nperms
           set.perms[set.keep]<-nperms
         }
         dimnames(gene.pvals)[[2]]<-gsub("stat", "perm.p", dimnames(gene.res)[[2]])
         dimnames(enrich.pvals)[[2]]<-gsub("stat", "perm.p", dimnames(enrich.res)[[2]])
 
         generes.tab<-cbind.data.frame(probeid=probes,
                                       gene.res,
                                       gene.pvals,
                                       nperms=num.perms)
         
         enrich.tab<-cbind.data.frame(setid=uset,
                                      enrich.res,
                                      enrich.pvals,
                                      nperms=set.perms)
         
         message(paste("Finished at: ",date()))
         
         res<-list(generes=generes.tab, setres=enrich.tab)
       }
       
       
       if (is.null(geneSet)) 
       {
         
         # Now compute results for observed data
         message(paste("Computing Stats for Observed Data:",date()))
         message(paste("Computing gene-specific statistics:",date()))
         if (is.null(strat.var)) strat<-NULL  
         else  strat<-phtype[,strat.var] 
         gene.res<-promise.genestat(Y,phtype,ph.pattern, strat, proj0=proj0)
         
         m.nph<-dim(gene.res)
         gene.pvals<-matrix(0,m.nph[1],m.nph[2])
         
         #set seed for sampling subjects
         message(paste("Set seed: seed = ",seed))
         set.seed(seed)
         
         probe.keep<-rep(T,nprobes)
         probe.done<-rep(F,nprobes)
         probe.ntail<-matrix(0,m.nph[1],m.nph[2])
         probe.apt.ntail<-rep(0,m.nph[1])
         num.perms<-rep(NA,nprobes)
         
         i<-0
         while((i<nperms)&&(any(probe.keep)))
         {
           i<-i+1
           if (is.null(strat.var)) perm.ind<-sample(nsamp.inc,replace=FALSE)
           else
           {
             strat<-X[,strat.var]
             ustrat<-sort(unique(strat))
             nstrat<-length(ustrat)
             ind<-length(strat)
             for (j in 1:nstrat)
             {
               this.strat<-(strat == ustrat[j])
               perm.ind[this.strat]<-sample((1:ind)[this.strat],replace=FALSE)
             }
           }
           perm.phtype<-phtype[perm.ind,]
           
           if (is.null(strat.var))  perm.strat<-NULL 
           else  perm.strat<- perm.phtype[,strat.var]
           
           gene.temp<-promise.genestat(Y[probe.keep, ], perm.phtype, ph.pattern, perm.strat, proj0=proj0)         
           probe.ntail[probe.keep,]<-probe.ntail[probe.keep,]+(abs(gene.temp)>=abs(gene.res[probe.keep,]))
           probe.apt.ntail[probe.keep]<- probe.ntail[probe.keep, m.nph[2]]
           
           probe.done[probe.keep]<-probe.apt.ntail[probe.keep]==max.ntail         
           num.perms[probe.done&probe.keep]<-i
           gene.pvals[probe.done&probe.keep,]<-probe.ntail[probe.done&probe.keep,]/i
           probe.keep[probe.done]<-F
         }
         
         if (any(probe.keep))
         {
           gene.pvals[probe.keep,]<-probe.ntail[probe.keep,]/nperms
           num.perms[probe.keep]<-nperms
         }
         dimnames(gene.pvals)[[2]]<-gsub("stat", "perm.p", dimnames(gene.res)[[2]])
         
         generes.tab<-cbind.data.frame(probeid=probes,
                                       gene.res,
                                       gene.pvals,
                                       nperms=num.perms)
         
         message(paste("Finished at: ",date()))
         
         res<-list(generes=generes.tab, setres=NULL)
       }
     } # end of nbperm 
     
     if (!nbperm)
     {
       if (!is.null(geneSet)) 
       {
         
         #Now extract probe and gene set names
         gs.probes<-geneset.data[,1]
         gs.sets<-geneset.data[,2]
         
         uset<-unique(gs.sets)
         nset<-length(uset)
         
         # Now compute results for observed data
         message(paste("Computing Stats for Observed Data:",date()))
         message(paste("Computing gene-specific statistics:",date()))
         if (is.null(strat.var)) strat<-NULL  
         else  strat<-phtype[,strat.var] 
         gene.res<-promise.genestat(Y,phtype,ph.pattern, strat, proj0=proj0)
         
         m.nph<-dim(gene.res)
         gene.pvals<-matrix(0,m.nph[1],m.nph[2])
         
         message(paste("Computing gene-set enrichment statistics:",date()))
         enrich.res<-avg.abs.genestat(gene.res,probes, geneset.data)
         
         ngs.nph<-dim(enrich.res)
         enrich.pvals<-matrix(0,ngs.nph[1],ngs.nph[2])
         
         #set seed for sampling subjects
         message(paste("Set seed: seed = ",seed))
         set.seed(seed)
         
         perms<-matrix(NA,nperms,nsamp.inc)
         if (is.null(strat.var)) for (i in 1:nperms) perms[i,]<-sample(nsamp.inc,replace=FALSE)
         else
         {
           strat<-endpt.data[,strat.var]
           ustrat<-sort(unique(strat))
           nstrat<-length(ustrat)
           ind<-length(strat)
           for (j in 1:nstrat)
           {
             this.strat<-(strat == ustrat[j])
             for (i in 1:nperms) perms[i,this.strat]<-sample((1:ind)[this.strat],replace=FALSE)
           }
         }
         
         message(paste("Computing Stats for Permuted Data: ",date()))    
         for (i in 1:nperms)
         {
           perm.ind<-unlist(perms[i,])
           perm.phtype<-phtype[perm.ind,]
           if (is.null(strat.var))  perm.strat<-NULL 
           else  perm.strat<- perm.phtype[,strat.var]
           
           gene.temp<-promise.genestat(Y,perm.phtype,ph.pattern, perm.strat, proj0=proj0)
           enrich.temp<-avg.abs.genestat(gene.temp,probes,geneset.data)
           
           enrich.pvals<-enrich.pvals+(abs(enrich.temp)>=abs(enrich.res))
           gene.pvals<-gene.pvals+(abs(gene.temp)>=abs(gene.res))
         }
         
         message(paste("Generating Result List: ",date()))
         gene.pvals<-gene.pvals/nperms
         enrich.pvals<-enrich.pvals/nperms

         dimnames(gene.pvals)[[2]]<-gsub("stat", "perm.p", dimnames(gene.res)[[2]])
         dimnames(enrich.pvals)[[2]]<-gsub("stat", "perm.p", dimnames(enrich.res)[[2]])
         
         generes.tab<-cbind.data.frame(probeid=probes,
                                       gene.res,
                                       gene.pvals)
         
         enrich.tab<-cbind.data.frame(setid=uset,
                                      enrich.res,
                                      enrich.pvals)
         
         message(paste("Finished at: ",date()))
         
         res<-list(generes=generes.tab, setres=enrich.tab)
       }
       
       if (is.null(geneSet))
       {
         # Now compute results for observed data
         message(paste("Computing Stats for Observed Data:",date()))
         message(paste("Computing gene-specific statistics:",date()))
         if (is.null(strat.var)) strat<-NULL  
         else  strat<-phtype[,strat.var] 
         gene.res<-promise.genestat(Y,phtype,ph.pattern, strat, proj0=proj0)
         
         m.nph<-dim(gene.res)
         gene.pvals<-matrix(0,m.nph[1],m.nph[2])
         
         #set seed for sampling subjects
         message(paste("Set seed: seed = ",seed))
         set.seed(seed)
         
         perms<-matrix(NA,nperms,nsamp.inc)
         if (is.null(strat.var)) for (i in 1:nperms) perms[i,]<-sample(nsamp.inc,replace=FALSE)
         else
         {
           strat<-endpt.data[,strat.var]
           ustrat<-sort(unique(strat))
           nstrat<-length(ustrat)
           ind<-length(strat)
           for (j in 1:nstrat)
           {
             this.strat<-(strat == ustrat[j])
             for (i in 1:nperms) perms[i,this.strat]<-sample((1:ind)[this.strat],replace=FALSE)
           }
         }
         
         message(paste("Computing Stats for Permuted Data: ",date()))    
         for (i in 1:nperms)
         {
           perm.ind<-unlist(perms[i,])
           perm.phtype<-phtype[perm.ind,]
           if (is.null(strat.var))  perm.strat<-NULL 
           else  perm.strat<- perm.phtype[,strat.var]
           
           gene.temp<-promise.genestat(Y,perm.phtype,ph.pattern, perm.strat, proj0=proj0)
           gene.pvals<-gene.pvals+(abs(gene.temp)>=abs(gene.res))
         }
         
         message(paste("Generating Result List: ",date()))
         
         gene.pvals<-gene.pvals/nperms
         dimnames(gene.pvals)[[2]]<-gsub("stat", "perm.p", dimnames(gene.res)[[2]])
         
         generes.tab<-cbind.data.frame(probeid=probes,
                                       gene.res,
                                       gene.pvals)
         message(paste("Finished at: ",date()))     
         res<-list(generes=generes.tab, setres=NULL)
       }       
       
     }

     return(res)
}
