#dir('../OscopeData/dataAdjacencyMatrixe02/',pattern='.csv')

library(igraph)
####GRAPH filename
genedet='../OscopeData/GenesLabels/Mapping/Water45.txt'
netdirList=c('../OscopeData/dataAdjacencyMatrixe02/','../OscopeData/dataAdjacencyMatrixePiDiv9/','../OscopeData/FirstData/')
outdirList=c('../Communities/Networks02/','../Communities/NetworksPiDiv9/','../Communities/FirstData/')
for (h in c(1:3)){
#netdir='../OscopeData/dataAdjacencyMatrixe02/'
#netdir='../OscopeData/dataAdjacencyMatrixePiDiv9/'
#netdir='../OscopeData/FirstData/'

#outdir='../Communities/Networks02/'
#outdir='../Communities/NetworksPiDiv9/'
#outdir='../Communities/FirstData/'

filename=paste(netdirList[h],'WaterfallbootstrapSequential_CPUDPSF_TF_2000Boot_data_s45res_C28G8442_edgeNetwork.csv',sep="")
print(filename)

data=read.csv(file=filename, header= TRUE,row.names=1)
G=graph_from_data_frame(data,directed=FALSE)
ww=get.edge.attribute(G,'cost')

###############WALKTRAP COMMUNITY 
wt=walktrap.community(G)
wwt=walktrap.community(G,weights=ww)
modularity(wwt)
#####modularity 0.3929114
modularity(wt)
#####modularity 0.5716109



#fastgreedy.community
wfg=fastgreedy.community(G,weights=ww)
fg=fastgreedy.community(G)
modularity(fg)
modularity(wfg)
# modularity(fg) 0.5980459
# modularity(wfg) 0.4243573

######################
###infomap
wim=infomap.community(G,v.weights=ww)
im=infomap.community(G)
modularity(im)
modularity(wim)
#modularity(im) 0.6025636
#modularity(wim) 0.6019861


######Label propagation
lp=label.propagation.community(G)
wlp=label.propagation.community(G,weights=ww)
modularity(lp)
modularity(wlp)
#modularity(lp) 0.5810269
#modularity(wlp) 0.3322666


compare(membership(fg), membership(wt),method="nmi")
compare(membership(wt), membership(im),method="nmi")

compare(membership(fg), membership(im),method="nmi")
compare(membership(fg), membership(wfg),method="nmi")

compare(membership(lp), membership(wlp),method="nmi")
compare(membership(lp), membership(im),method="nmi")
compare(membership(lp), membership(wt),method="nmi")
compare(membership(lp), membership(fg),method="nmi")
####all the scores between 0.96 and 0.97, maximum nmi is 1! so all the methods are in a quite high agreement#####
 max(sizes(wt))
#[1] 739
max(sizes(im))
#[1] 194
 max(sizes(fg))
#[1] 492
 max(sizes(lp))
#[1] 1985



Mwt=as.matrix(membership(wt))
Mim=as.matrix(membership(im))
Mfg=as.matrix(membership(fg))



genes=read.table(genedet,sep='\t',header=TRUE, quote="")
merMwt=merge(Mwt,genes,by.x=0,by.y='mgi_symbol',all.x=TRUE)
merMim=merge(Mim,genes,by.x=0,by.y='mgi_symbol',all.x=TRUE)
merMfg=merge(Mfg,genes,by.x=0,by.y='mgi_symbol',all.x=TRUE)
cn=c('mgi_symbol','Community Id')
colnames(merMwt)[1:2]=cn
colnames(merMim)[1:2]=cn
colnames(merMfg)[1:2]=cn

write.table(merMwt, file=paste(outdirList[h],'W45walktrap.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
write.table(merMim, file=paste(outdirList[h],'W45infomap.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
write.table(merMfg, file=paste(outdirList[h],'W45fastgreedy.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)


}