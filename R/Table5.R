#/Users/luisacutillo/Manchester/Alexis/Oscope++/OscopeRcode/casestudy/H1
library('igraph')
library('biomaRt')
library('neat')
library('Oscope')
mypath='~/Manchester/Alexis/Oscope++/'
source(paste(mypath,'OscopeRcode/CustomLibrary.R',sep=''))

#mypath='/home/luisa/Manchester/Alexis/Oscope++/'
mypath='/Users/luisacutillo/Manchester/Alexis/Oscope++/'
customlib=paste(mypath,'OscopeRcode/CustomLibrary.R',sep='')

source(customlib)

workingdir=paste(mypath,'OscopeRcode/casestudy/H1/results/',sep='')
outDir='Communities/'

filelist=dir(path=workingdir,pattern='mergedH1_N1000_g18_TFTrue_SummaryPartition')
beforefilt=file.path('~/Manchester/Alexis/Oscope++/OscopeRcode/casestudy/H1/mergedH1.csv')
beforedata=read.csv(file=beforefilt, header= TRUE,row.names=1)
originalGenes=row.names(beforedata)

#####for loop


filea=paste(workingdir,filelist[1],sep='')
ppp=strsplit(filelist[1],split=".csv")
basefilea=ppp[[1]]
A=read.csv(file=filea,header=TRUE,sep=',')

G=graph_from_data_frame(A[,c(2,3,4)],directed =FALSE)

# #####apply the KM Oscope Method
# MMM=get.adjacency(G,attr = 'cost')
# ##eps2=get.edge.attribute(G,'cost')
# ##min(eps2) ##this is greather then zero
# cut=-log10(min(MMM[which(as.matrix(MMM)!=0)]))
# MMM[which(as.matrix(MMM)==0)]=2*max(MMM)
# SimiMat=-log10(MMM)
# MyRes=list(SimiMat=as.matrix(SimiMat))
# KMmyRes <- OscopeKM(MyRes, cut = cut, maxK = 100, minSize = 1, maxSize = 1000)
# print(KMmyRes)
# 
###############WALKTRAP COMMUNITY 
wt=walktrap.community(G,steps=4)
modularity(wt)
mwt=membership(wt)
sz=sizes(wt)
Mwt=as.matrix(membership(wt))
genes=row.names(Mwt)

caseE='ensembl_gene_id' #GENESYMBOL'
caseS="hgnc_symbol"###in case of mouse mgi_symbol
case=caseS
#speces="mmusculus_gene_ensembl"
speces="hsapiens_gene_ensembl"

#mart <- useDataset(speces, useMart("ensembl"))
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
mart = useMart("ensembl")
mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")




CC=getBM(filters="go",
      values=c('GO:0007049'), mart=mart,
      attributes= c(caseE,caseS,"description"))
G_list <- getBM(filters= case, attributes= c(caseE,
                                            caseS,"description"),values=genes,mart= mart)
commonCCbefore=intersect(CC$hgnc_symbol,originalGenes)

commonCC=intersect(CC$hgnc_symbol,V(G)$name)


details=as.matrix(G_list)
colnames(Mwt)='CommunityID'
merMwt=Mwt#merge(details,Mwt,by.x=case,by.y=0,all.y=TRUE)
write.table(merMwt, file=paste(workingdir,outDir,basefilea,'Comm.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
####merge communities with Cell Cycle related genes
CC=unique(CC$hgnc_symbol)
CC=as.data.frame(CC)
names(CC)='hgnc_symbol'
merMwt=data.frame(hgnc_symbol=row.names(merMwt),merMwt)

CCmerMwt=merge(CC,merMwt,by.x=caseS,by.y=caseS)
write.table(CCmerMwt, file=paste(workingdir,outDir,basefilea,'CC_Comm.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)

vecComCC=unique(CCmerMwt$CommunityID)
sizeCCcom=vecComCC
EnrichmentTest=vecComCC
GlobalRatio=vecComCC
RelativeRatio=vecComCC

white=length(unique(CCmerMwt$hgnc_symbol))#TotalCommonTargets
black=length(unique(merMwt$hgnc_symbol))-white
validation1=NULL
validation2=NULL
for (h in c(1:length(vecComCC))){
  
  ####neat testing for edges
  genlist1=unique(CCmerMwt$hgnc_symbol[which(CCmerMwt$CommunityID==vecComCC[h])])
  #myList[[length(myList)+1]] <- list(genlist1)
  genlist2=unique(merMwt$hgnc_symbol[which(merMwt$CommunityID==vecComCC[h])])
  alist = list('set 1' = genlist1)
  names(alist)=vecComCC[h]
  blist= list('set 2'=genlist2)
  test1 = neat(alist = alist, blist = alist, network = G,
               nettype = 'undirected', nodes = V(G)$name, alpha = 0.1)
  test2 = neat(alist = alist, blist = blist, network = G,
               nettype = 'undirected', nodes = V(G)$name, alpha = 0.1)
  validation1=rbind(validation1,print(test1))
  validation2=rbind(validation2,print(test2))
  
  ####hypergeometric testing for vertex
  sizeCCcom[h]=length(which(CCmerMwt$CommunityID==vecComCC[h]))
  x=sizeCCcom[h]#intersection
  k=length(which(merMwt$CommunityID==vecComCC[h]))#commsize
  EnrichmentTest[h]=1-phyper(x,white,black,k)
  GlobalRatio[h]=x/white
  RelativeRatio[h]=x/k
}
BHE=p.adjust(EnrichmentTest,method='BH')
outCC=cbind(vecComCC,sizeCCcom,GlobalRatio,RelativeRatio,EnrichmentTest,BHE)
colnames(outCC)=c('CC enriched Comm ID', 'size CC Overlapping','Global CC Ratio','Relative CC Ratio','CC Overlapping Significance','CC sig BH corr')
CCvalidation=validation1



IdGE10=which(sz>1)
n10=length(IdGE10)
MatSummary=matrix('NA',n10,6)
colnames(MatSummary)=c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max." )
NumberOfEdges=matrix('NA',n10,1)
CommDensity=matrix('NA',n10,1)
RelativeCommDensity=matrix('NA',n10,1)
PsigCom=matrix('NA',n10,1)
TotalEdges=ecount(G)
for (i in c(1:n10)){
  actCom=IdGE10[i]
  actNodes=which(Mwt==actCom)
  actNumVertex=length(actNodes)
  subG=induced_subgraph(G,actNodes)
  actEdge=E(subG)
  actNumEdges=ecount(subG)
  NumberOfEdges[i]=actNumEdges
  CommDensity[i]=actNumEdges/TotalEdges
  RelativeCommDensity[i]=2*actNumEdges/(actNumVertex^2-actNumVertex)
  
  actCost=get.edge.attribute(subG,'cost',index=actEdge)
  vsumAct=as.vector(summary(actCost))
  MatSummary[i,]=vsumAct
  res=community.significance.test(G, actNodes)
  PsigCom[i]=res$p.value
}
sizes10=sz[IdGE10]
sizes10=data.frame(Community=names(sizes10),NumberOfNodes= as.vector(sizes10))

Report10=data.frame(sizes10,NumberOfEdges=NumberOfEdges,Significance=PsigCom,Density=CommDensity,RelativeDensity=RelativeCommDensity,MatSummary)
Report10sig=Report10[which(as.numeric(as.vector(Report10$Significance))<0.01),]
write.table(Report10sig, file=paste(workingdir,'Report1sig','.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)

basicInfo=data.frame(sizes10,NumberOfEdges=NumberOfEdges,Significance=PsigCom)
basicInfo=basicInfo[which(as.numeric(as.vector(basicInfo$Significance))<0.09),]

valCC=merge(outCC,basicInfo,by.y='Community',by.x='CC enriched Comm ID')

valCC=merge(valCC,CCvalidation,by.x='CC enriched Comm ID',by.y='A',all.x=TRUE)
write.table(valCC, file=paste(workingdir,'Report1CC','.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)

indSigCC=which(as.numeric(as.vector(valCC$Significance))<0.09)
write.table(valCC[indSigCC,], file=paste(workingdir,'Neat/',basefilea,'GlobalNeatCC','.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)

valCC=valCC[indSigCC,]
valCC=valCC[valCC$pvalue<0.05,]
idFinalCC=valCC$`CC enriched Comm ID`[valCC$`size CC Overlapping`>1]
outNames <- paste0("cluster", idFinalCC)
CCmerMwt$hgnc_symbol[CCmerMwt$CommunityID==idFinalCC[1]]
out=vector("list",length(idFinalCC))
names(out)=outNames
for (jj in c(1:length(idFinalCC))){
  out[[jj]]=CCmerMwt$hgnc_symbol[CCmerMwt$CommunityID==idFinalCC[jj]]
}

ns=dim(Report10sig)[1]
out10Sig=vector("list",ns)
outNamesSig <- paste0("cluster", Report10sig$Community)
names(out10Sig)=outNamesSig
for (jj in c(1:ns)){
  out10Sig[[jj]]=merMwt[merMwt$CommunityID==Report10sig$Community[jj],]
}
#######

save(out10Sig,Report10,Report10sig,out,CCmerMwt,valCC,CC,file=file.path('~/Manchester/Alexis/Oscope++/OscopeRcode/casestudy/H1/','ClusterH1merged.Rdata'))
#preparing input for cytoscape

#######
#preparing input for cytoscape
outnet='/Users/luisacutillo/Manchester/Alexis/Oscope++/OscopeRcode/Cytoscape/Hsubnet.txt'
outnetN='/Users/luisacutillo/Manchester/Alexis/Oscope++/OscopeRcode/Cytoscape/HsubnetNames.txt'

outlabel='/Users/luisacutillo/Manchester/Alexis/Oscope++/OscopeRcode/Cytoscape/Hlabels.txt'
outlabelN='/Users/luisacutillo/Manchester/Alexis/Oscope++/OscopeRcode/Cytoscape/HlabelsNames.txt'

ids=as.numeric(as.vector(basicInfo$Community))
vetnodes=NULL
vetcom=NULL
vetNames=NULL
namesV=V(G)
for(i in c(1:length(ids))){
  dum=which(Mwt==ids[i])
  vetnodes=c(vetnodes,dum)
  vetNames=c(vetNames,namesV[dum])
  vetcom=c(vetcom,rep(ids[i],length(dum)))
}
subgN<-induced.subgraph(G,names(vetNames))
subg1<-induced.subgraph(G,vetnodes)
NodesCom=data.frame(NODESid=vetnodes,COM=vetcom)
NodesComNames=data.frame(NodesName=names(vetNames),NODESid=vetnodes,COM=vetcom)
write_graph(subg1, outnet, "edgelist")
write_graph(subgN, outnetN, "ncol")
write.table(NodesComNames,file=outlabelN,sep='\t')
write.table(NodesCom,file=outlabel,sep='\t')
# 
# Report10=merge(Report10,outCC,by.x='Community',by.y='CC enriched Comm ID',all.x=TRUE)
# #outmir294TmerMwt
# #outlet7TmerMwt
# indSig=which(as.numeric(as.vector(Report10$Significance))<0.09)
# 
# 
# 
# 
# 
# write.table(Report10[indSig,], file=paste(workingdir,outDir,basefilea,'Report5_Comm.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)


# ############plot only the significant communities
# comm=wt
# Net=G
# V(Net)$comm=membership(comm)
# subG=induced_subgraph(Net,membership(comm)%in%Report10sig$Community)
# Net=subG
# 
# V(Net)$size = 3
# nc=length(Report10sig$Community)
# newvals=c(1:length(Report10sig$Community))
# V(Net)$color=newvals[as.factor(V(Net)$comm)]
# pdf('~/Manchester/Alexis/Oscope++/OscopeRcode/casestudy/H1/H1SigClusNet.pdf')
# plot(Net, layout=layout.fruchterman.reingold, vertex.color=V(Net)$color, vertex.size = V(Net)$size, vertex.label = NA, edge.arrow.size = 0.5)
# dev.off()
# 
# #######contract comunities
# 
# Net=subG
# 
# V(Net)$size = 1
# E(Net)$count = 1
# nc=length(Report10sig$Community)
# newvals=c(1:nc)
# V(Net)$comm=newvals[as.factor(V(Net)$comm)]
# comm.graph <- contract.vertices(Net, V(Net)$comm, vertex.attr.comb=list(size="sum", "ignore"))
# comm.graph <- simplify(comm.graph, remove.loops=TRUE, edge.attr.comb=list(count="sum", "ignore"))
# M=max(E(comm.graph)$count)
# 
# vectsizes=V(comm.graph)$size
# scale=50
# V(comm.graph)$name <- LETTERS[1:vcount(comm.graph)]#as.character(Report10sig$Community)
# ###############Set the best colors range you prefer
# V(comm.graph)$color <- c(1:vcount(comm.graph))
# V(comm.graph)$size <- vectsizes/max(vectsizes)*scale
# E(comm.graph)$weight <-E(comm.graph)$count
# E(comm.graph)$length <-M-E(comm.graph)$count#1/E(comm.graph)$count
# pdf('~/Manchester/Alexis/Oscope++/OscopeRcode/casestudy/H1/H1ContractedSigNet.pdf')
# 
# plot(comm.graph,layout=layout.fruchterman.reingold,vertex.label.dist=0)
# 
# dev.off()
# 
# #colrs <- c("gray50", "tomato", "gold", "blue", "yellow")
# 
# #V(g)$color <- colrs[V(g)$name]
# ###plot the edges width according to the number of edges etween communities
# #
# 
