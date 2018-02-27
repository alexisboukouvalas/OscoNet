library('igraph')
library('biomaRt')
library('neat')
library('Oscope')
mypath='./'

customlib=paste(mypath,'CustomLibrary.R',sep='')

source(customlib)

workingdir= mypath#paste(mypath,'OscopeRcode/casestudy/H1/results/',sep='')
#outDir=mypath#'Communities/'

filea=file.path('./toy.csv')



A=read.csv(file=filea,header=TRUE,sep=',')

G=graph_from_data_frame(A[,c(2,3,4)],directed =FALSE)

npresent=length(get.vertex.attribute(G,"name"))
###############WALKTRAP COMMUNITY 
wt=walktrap.community(G,steps=4)
modularity(wt)
mwt=membership(wt)
sz=sizes(wt)
Mwt=as.matrix(membership(wt))
genes=row.names(Mwt)

colnames(Mwt)='CommunityID'
merMwt=Mwt#merge(details,Mwt,by.x=case,by.y=0,all.y=TRUE)
write.table(merMwt, file=paste(workingdir,filea,'Comm.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)

##we are excluding singletons and performinf a significance test on the communities
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

#This contains all the basic info about the communities with at least 2 nodes
#basicInfo=data.frame(sizes10,NumberOfEdges=NumberOfEdges,Significance=PsigCom)
#basicInfo=basicInfo[which(as.numeric(as.vector(basicInfo$Significance))<0.09),]
