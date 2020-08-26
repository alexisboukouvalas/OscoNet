source('GOEA.R')

filecomm='../Communities/Networks02/W45fastgreedy.txt'
thresh=80

AAA=read.delim(filecomm)
uc=unique(AAA$Community.Id)
nuc=length(uc)
universe=as.vector(AAA$mgi_symbol)
results=NULL
for (i in c(1:nuc)){
	dum=which(AAA$Community.Id%in%uc[i])
	ld=length(dum)
	if(ld>=thresh){
		print(i)
		group= universe[dum]
		gene_sets=list(group)
		names(gene_sets)=paste("comm",i,sep='')
		enrichLists <- do.GOEA(gene_sets,universe)
		results=list(results,get.df(enrichLists,th = 2,th_fdr = 2))
		
	}
	
	
	
}