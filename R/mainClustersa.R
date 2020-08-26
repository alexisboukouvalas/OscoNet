#dir('../OscopeData/dataAdjacencyMatrixe02/',pattern='.csv')

library(igraph)
####GRAPH filename
genedet='../OscopeData/GenesLabels/Mapping/SVZa.txt'
netdirList=c('../OscopeData/dataAdjacencyMatrixe02/','../OscopeData/dataAdjacencyMatrixePiDiv9/','../OscopeData/FirstData/')
outdirList=c('../Communities/Networks02/','../Communities/NetworksPiDiv9/','../Communities/FirstData/')

#netdir='../OscopeData/dataAdjacencyMatrixe02/'
#netdir='../OscopeData/dataAdjacencyMatrixePiDiv9/'
#netdir='../OscopeData/FirstData/'

#outdir='../Communities/Networks02/'
#outdir='../Communities/NetworksPiDiv9/'
#outdir='../Communities/FirstData/'

for(h in c(1:3)){
	
	
filename=paste(netdirList[h],'SVZbootstrap_CPUDPSF_TF_2000Boot_data_aNSC_C32G7883_edgeNetwork.csv',sep="")
print(filename)
data=read.csv(file=filename, header= TRUE,row.names=1)
G=graph_from_data_frame(data,directed=FALSE)
ww=get.edge.attribute(G,'cost')

###############Leading Eignevectors
le=cluster_leading_eigen(G)
modularity(le) #0.658891
wle=cluster_leading_eigen(G, weights=ww)
###############WALKTRAP COMMUNITY 
wt=walktrap.community(G)
modularity(wt)
mwt=membership(wt)
#####modularity  0.657515

# ########plot the sizes of the communities
# hist(sizes(wt),300)
# #####as you can see there are few big ones and many single nodes.

# m=10
# x=which(sizes(wt)>=m)
# ind=NULL
# clust=NULL
# for (i in c(1:length(x))){
	# pos=which(mwt== x[i])
	# ind=c(ind,pos )
	# clust=c(clust,rep(i,1,length(pos)))
# }
# #clust=meb[ind]

# subg <- induced.subgraph(G, ind)

# colors <- rainbow(max(clust))
# #plot(subg,vertex.color=colors[clust], 
 # #    layout=layout.fruchterman.reingold, vertex.size=3, vertex.label=NA)

###############try fast greedy modularity optimization method

#fastgreedy.community
wfg=fastgreedy.community(G,weights=ww)
fg=fastgreedy.community(G)
modularity(fg)#0.684239

modularity(wfg)#0.6972934

######################
###infomap
wim=infomap.community(G,v.weights=ww)
im=infomap.community(G)
modularity(im)#0.6784341
modularity(wim)#0.6784341

######Label propagation
#lp=label.propagation.community(G)
#wlp=label.propagation.community(G,weights=ww)
#modularity(lp)#0.5590939
#modularity(wlp)#0.5873755

# compare(membership(fg), membership(wt),method="nmi")
# compare(membership(wt), membership(im),method="nmi")

# compare(membership(fg), membership(im),method="nmi")
# compare(membership(lp), membership(im),method="nmi")
# compare(membership(lp), membership(wt),method="nmi")
# compare(membership(lp), membership(fg),method="nmi")

#%%%all very similar but wt and fg are the closest 0.95

Mwt=as.matrix(membership(wt))
Mim=as.matrix(membership(im))
Mfg=as.matrix(membership(fg))


genes=read.table(genedet,sep="\t",header=TRUE,quote="")
merMwt=merge(Mwt,genes,by.x=0,by.y='ensembl_gene_id',all.x=TRUE)
merMim=merge(Mim,genes,by.x=0,by.y='ensembl_gene_id',all.x=TRUE)
merMfg=merge(Mfg,genes,by.x=0,by.y='ensembl_gene_id',all.x=TRUE)
cn=c('ensembl_gene_id','Community Id')
colnames(merMwt)[1:2]=cn
colnames(merMim)[1:2]=cn
colnames(merMfg)[1:2]=cn

write.table(merMwt, file=paste(outdirList[h],'Sawalktrap.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
write.table(merMim, file=paste(outdirList[h],'Sainfomap.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
write.table(merMfg, file=paste(outdirList[h],'Safastgreedy.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)

}
