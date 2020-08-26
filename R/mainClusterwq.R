library(igraph)

####GRAPH filename
genedet='../OscopeData/GenesLabels/Mapping/Water123.txt'
netdirList=c('../OscopeData/dataAdjacencyMatrixe02/','../OscopeData/dataAdjacencyMatrixePiDiv9/','../OscopeData/FirstData/')
outdirList=c('../Communities/Networks02/','../Communities/NetworksPiDiv9/','../Communities/FirstData/')
for (h in c(1:3)){
filename=paste(netdirList[h],'WaterfallbootstrapSequential_CPUDPSF_TF_2000Boot_data_s123res_C73G9200_edgeNetwork.csv',sep="")

print(filename)

#netdir='../OscopeData/dataAdjacencyMatrixe02/'
#netdir='../OscopeData/dataAdjacencyMatrixePiDiv9/'
#netdir='../OscopeData/FirstData/'


#outdir='../Communities/Networks02/'
#outdir='../Communities/NetworksPiDiv9/'
#outdir='../Communities/FirstData/'






data=read.csv(file=filename, header= TRUE,row.names=1)
G=graph_from_data_frame(data,directed=FALSE)
ww=get.edge.attribute(G,'cost')


###############WALKTRAP COMMUNITY 
wt=walktrap.community(G)
wwt=walktrap.community(G,weights=ww)
modularity(wwt)
modularity(wt)

#modularity(wwt) 0.3090623
#modularity(wt) 0.3156759

#####before modularity 0.3963771



#fastgreedy.community
wfg=fastgreedy.community(G,weights=ww)
fg=fastgreedy.community(G)
modularity(fg)
modularity(wfg)

#modularity(fg) 0.3909607
#modularity(wfg) 0.3855991


######################
###infomap
wim=infomap.community(G,v.weights=ww)
im=infomap.community(G)
modularity(im)
modularity(wim)

#modularity(im) 0.3747024
#modularity(wim) 0.3747024

######Label propagation
lp=label.propagation.community(G)
wlp=label.propagation.community(G,weights=ww)
modularity(lp)
modularity(wlp)

#modularity(lp) 0.1322383
#modularity(wlp) 0.1283502

compare(membership(fg), membership(wt))
compare(membership(wt), membership(im))

compare(membership(fg), membership(im))
compare(membership(fg), membership(wfg))

compare(membership(lp), membership(wlp))
compare(membership(lp), membership(im))
compare(membership(lp), membership(wt))
compare(membership(lp), membership(fg))

####all the scores between 0.79 and 0.81, maximum nmi is 1! so all the methods are in a quite high agreement#####


#############################################All methods show a  lower modularity in w123 with respet to w45!
genes=read.table(genedet,sep='\t',header=TRUE, quote="")
Mwt=as.matrix(membership(wt))
Mim=as.matrix(membership(im))
Mfg=as.matrix(membership(fg))


merMwt=merge(Mwt,genes,by.x=0,by.y='mgi_symbol',all.x=TRUE)
merMim=merge(Mim,genes,by.x=0,by.y='mgi_symbol',all.x=TRUE)
merMfg=merge(Mfg,genes,by.x=0,by.y='mgi_symbol',all.x=TRUE)
cn=c('mgi_symbol','Community Id')
colnames(merMwt)[1:2]=cn
colnames(merMim)[1:2]=cn
colnames(merMfg)[1:2]=cn

 
write.table(merMwt, file=paste(outdirList[h],'W123walktrap.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
write.table(merMim, file=paste(outdirList[h],'W123infomap.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
write.table(merMfg, file=paste(outdirList[h],'W123fastgreedy.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)

}