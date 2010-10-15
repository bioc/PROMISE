`PROMISE` <-
function(exprSet, geneSet=NULL, promise.pattern, strat.var=NULL, seed=13, nperms=10000)

{
     # extract data components from objects
     ph.pattern<-promise.pattern
     endpt.data<-pData(phenoData(exprSet))
     array.data<-exprs(exprSet)
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

     if (is.vector(phtype)) phtype<-matrix(phtype,length(phtype),1)
     
     # Now extract probe set names
     probes<-dimnames(array.data)[[1]]
     nprobes<-length(probes)
     
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
         gene.res<-promise.genestat(Y,phtype,ph.pattern, strat)

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
                   
             gene.temp<-promise.genestat(Y,perm.phtype,ph.pattern, perm.strat)
             enrich.temp<-avg.abs.genestat(gene.temp,probes,geneset.data)

             enrich.pvals<-enrich.pvals+(abs(enrich.temp)>=abs(enrich.res))
             gene.pvals<-gene.pvals+(abs(gene.temp)>=abs(gene.res))
         }
     
         message(paste("Generating Result List: ",date()))
         gene.pvals<-gene.pvals/nperms
         enrich.pvals<-enrich.pvals/nperms

         generes.tab<-cbind.data.frame(probeid=probes,
                                   gene.res,
                                   perm.p=gene.pvals)

         enrich.tab<-cbind.data.frame(setid=uset,
                                  enrich.res,
                                  perm.p=enrich.pvals)

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
        gene.res<-promise.genestat(Y,phtype,ph.pattern, strat)

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
                   
            gene.temp<-promise.genestat(Y,perm.phtype,ph.pattern, perm.strat)
            gene.pvals<-gene.pvals+(abs(gene.temp)>=abs(gene.res))
        }
     
        message(paste("Generating Result List: ",date()))

        gene.pvals<-gene.pvals/nperms
  
        generes.tab<-cbind.data.frame(probeid=probes,
                                   gene.res,
                                   perm.p=gene.pvals)
        message(paste("Finished at: ",date()))     
        res<-list(generes=generes.tab, setres=NULL)
     }
     return(res)
}

