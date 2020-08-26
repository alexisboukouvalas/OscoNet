###MainClusterWhitfild
library('igraph')
library('biomaRt')
library('neat')
library('Oscope')
load('~/Manchester/Alexis/Oscope++/OscopeRcode/casestudy/WhitfieldStandardOscope/StandardOscope.RData')

#mypath='/home/luisa/Manchester/Alexis/Oscope++/'
mypath='~/Manchester/Alexis/Oscope++/'
source(paste(mypath,'OscopeRcode/CustomLibrary.R',sep=''))

workingdir=paste(mypath,'OscopeRcode/casestudy/Whitfield/',sep='')
outDir=paste(mypath,'OscopeRcode/casestudy/Whitfield/','results/',sep='')



filename=file.path(workingdir,'filt0.9Data_N2000_g50_TFTrue_SummaryPartition.csv')
ppp=strsplit(filename,split=".csv")
basefilea=ppp[[1]]
#filename='filt0.9Data_N1000_g18_TFTrue_SummaryPartition.csv'
print(filename)
data=read.csv(file=filename, header= TRUE,row.names=1)
G=graph_from_data_frame(data,directed=FALSE)
#ww=get.edge.attribute(G,'cost')

###############WALKTRAP COMMUNITY 
wt=walktrap.community(G)
modularity(wt)
mwt=membership(wt)
sz=sizes(wt)
#sz[which(sz>=10)]
Mwt=as.matrix(membership(wt))

genes=read.delim(file.path('~/Manchester/Alexis/Oscope++/OscopeRcode/casestudy/Whitfield/ImageIdConversions.txt'),sep="\t",header=TRUE,quote="")
details=genes[,c(1:3)]
colnames(details)=c('CloneID','Ref','GeneSymbol')
colnames(Mwt)='CommunityID'

merMwt=merge(details,Mwt,by.x='CloneID',by.y=0,all.y=TRUE)
write.table(merMwt, file=file.path(outDir,'witfildCommG50.txt'),sep='\t',quote=FALSE,row.names=FALSE)
CC=read.delim(file.path('~/Manchester/Alexis/Oscope++/OscopeRcode/casestudy/Whitfield/CellCycleGeneList_1134.txt'))
#CCmerMwt=merge(CC,merMwt,by.x='CLONEID',by.y='CloneID')
CCmerMwt=merge(CC,Mwt,by.x='CLONEID',by.y=0)

write.table(CCmerMwt, file=file.path(outDir,'CCwitfildCommG50.txt'), sep='\t',quote=FALSE,row.names=FALSE)


vecComCC=unique(CCmerMwt$CommunityID)#community arricchite per CC
sizeCCcom=vecComCC
EnrichmentTest=vecComCC
GlobalRatio=vecComCC
RelativeRatio=vecComCC

white=length(unique(CCmerMwt$CLONEID))#TotalCommonTargets 360 CC presenti nella rete
black=length(unique(merMwt$CloneID))-white #Total Non CC, tutti gli altri
validation1=NULL
validation2=NULL

for (h in c(1:length(vecComCC))){
  
  ####neat testing for edges
  genlist1=CCmerMwt$CLONEID[which(CCmerMwt$CommunityID==vecComCC[h])]#CC genes in  cluster vecComCC[h]
  #myList[[length(myList)+1]] <- list(genlist1)
  genlist2=merMwt$CloneID[which(merMwt$CommunityID==vecComCC[h])]#All the genes in Cluster vecComCC[h]
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
outCC=cbind(outCC,CCvalidation)

##########################################
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

Report10=data.frame(sizes10,Significance=PsigCom,Density=CommDensity,RelativeDensity=RelativeCommDensity,MatSummary)
Report10sig=Report10[which(as.numeric(as.vector(Report10$Significance))<0.01),]
###########
###basicInfo contains the basic info for the Significant communities:
#ID, Significance, Number of nodes, Density and relative density
#we have 8 communities (25  38  61  70  74  92  124 139) significant
basicInfo=Report10sig[,c(1:5)]#data.frame(sizes10,NumberOfNodes=NumberOfNodes,NumberOfEdges=NumberOfEdges,Significance=PsigCom)
valCC=merge(outCC,basicInfo,by.y='Community',by.x=1)
AAA=merge(basicInfo,outCC,by.x='Community',by.y=1,all.x=TRUE)

Table1=data.frame(Size=as.vector(AAA$NumberOfNodes),CC=as.vector(AAA$'size CC Overlapping'),Density=as.vector(AAA$Density),Realative_density=as.vector(AAA$RelativeDensity),Hypergeometic=as.vector(AAA$'CC sig BH corr'), Significance=AAA$Significance, Neat=as.vector(AAA$adjusted_p))

Table1=Table1[order(Table1$Size),]