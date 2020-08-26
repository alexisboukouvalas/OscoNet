
library(GO.db)
######################
library(biomaRt)
library(limma)
Dirlist=c('../Communities/Networks02/','../Communities/NetworksPiDiv9/','../Communities/FirstData/')
#file='W45fastgreedy.txt'
file='Safastgreedy.txt'

thresh=8


for (h in c(1:3)){

Dircom= Dirlist[h]
#Dircom='../Communities/Networks02/'
filecomm=paste(Dircom,file, sep='')
outGo=paste(Dircom,'GO/Go_',file,sep='')
outKegg=paste(Dircom,'Kegg/Kegg_',file,sep='')




AAA=read.delim(filecomm)
uc=unique(AAA$Community.Id)
nuc=length(uc)
universe=as.vector(AAA$mgi_symbol)



mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))


 
 



G_list <- getBM(filters= "mgi_symbol", attributes= c("mgi_symbol",'entrezgene'),values=universe,mart= mart)
universeEntr=G_list$entrezgene

gene_sets=list()
names=NULL
sizes=NULL

for (i in c(1:nuc)){
	dum=which(AAA$Community.Id==uc[i])
	ld=length(dum)
	if(ld>=thresh){
		print(uc[i])
		print(ld)
		group= AAA$mgi_symbol[dum]
		a=merge(group,G_list,by.x=1,by.y=1)
		gene_sets[length(gene_sets)+1]=list(a$entrezgene)
		names=c(names,paste("comm",uc[i],sep=''))
		sizes=c(sizes,ld)
		}
		}
		names(gene_sets)=names
		ngs=length(names)
		
		GoRes=list()
		KeggRes=list()
		# for (i in c(1:ngs)){
			# enrichmentGO <- goana(gene_sets[[i]], universe=universeEntr,species='Mm')
# #enrichmentKEGG
			# GoRes[length(GoRes)+1]=list(enrichmentGO)
# #KeggRes[length(KeggRes)+1]=list(enrichmentKEGG)

		# }
        GoRes = goana(gene_sets,universe=universeEntr,species='Mm')
        KeggRes = kegga(gene_sets,universe=universeEntr,species='Mm')

write.table(cbind(rownames(GoRes),GoRes),row.names=FALSE,file=outGo,sep='\t',quote=FALSE)
write.table(cbind(rownames(KeggRes), KeggRes),file=outKegg,sep='\t',quote=FALSE,row.names=FALSE)
}