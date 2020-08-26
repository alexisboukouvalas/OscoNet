library(igraph)
####GRAPH filename
genedet='../OscopeData/GenesLabels/Mapping/SVZq.txt'
netdirList=c('../OscopeData/dataAdjacencyMatrixe02/','../OscopeData/dataAdjacencyMatrixePiDiv9/','../OscopeData/FirstData/')
outdirList=c('../Communities/Networks02/','../Communities/NetworksPiDiv9/','../Communities/FirstData/')

for (h in c(1:3)){



filename=paste(netdirList[h],'SVZbootstrap_CPUDPSF_TF_2000Boot_data_qNSC_C60G7883_edgeNetwork.csv',sep="")
print(filename)

data=read.csv(file=filename, header= TRUE,row.names=1)
G=graph_from_data_frame(data,directed=FALSE)
ww=get.edge.attribute(G,'cost')

###############Leading Eignevectors
le=cluster_leading_eigen(G)
modularity(le) #0.4244946
wle=cluster_leading_eigen(G, weights=ww)
modularity(wle)#0.439686
###############WALKTRAP COMMUNITY 
wt=walktrap.community(G)
wwt=walktrap.community(G,weights=ww)
modularity(wwt)#0.06753387
modularity(wt)#0.3997952
#####



#fastgreedy.community
wfg=fastgreedy.community(G,weights=ww)
fg=fastgreedy.community(G)
modularity(fg)#[1] 0.3679955

modularity(wfg)#0.4412479


######################
###infomap
wim=infomap.community(G,v.weights=ww)
im=infomap.community(G)
modularity(im)#0.3743114
modularity(wim)#0.4072138

######Label propagation
lp=label.propagation.community(G)
wlp=label.propagation.community(G,weights=ww)
modularity(lp)#0.2040412
modularity(wlp)# 0.1595253


compare(membership(fg), membership(wt),method="nmi")
compare(membership(wt), membership(im),method="nmi")
compare(membership(fg), membership(im),method="nmi")
compare(membership(lp), membership(im),method="nmi")
compare(membership(lp), membership(wt),method="nmi")
compare(membership(lp), membership(fg),method="nmi")

#%%%all  similar (0.75-0.85) but wt and im are the closest  0.8565346



max(sizes(wt))
#[1] 1710
max(sizes(im))
#[1] 2835
max(sizes(fg))
#[1] 2871
max(sizes(le))
#[1] 3122
max(sizes(lp))
#[1] 7225

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


write.table(merMwt, file=paste(outdirList[h],'Sqwalktrap.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
write.table(merMim, file=paste(outdirList[h],'Sqinfomap.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
write.table(merMfg, file=paste(outdirList[h],'Sqfastgreedy.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)

}
