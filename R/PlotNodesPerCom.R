
library(igraph)
####GRAPH filename
filename='../OscopeData/FirstData/WaterfallbootstrapSequential_CPUDPSF_TF_2000Boot_data_s45res_C28G8442_edgeNetwork.csv'
genedet='../OscopeData/GenesLabels/Mapping/Water45.txt'

data=read.csv(file=filename, header= TRUE,row.names=1)
G=graph_from_data_frame(data,directed=FALSE)
ww=get.edge.attribute(G,'cost')


###############WALKTRAP COMMUNITY 
wt=walktrap.community(G)
wwt=walktrap.community(G,weights=ww)
modularity(wwt)
modularity(wt)
 comm=wt
 Net=G
 V(Net)$comm=membership(comm)
 Thresh=20
 icomm=which(sizes(comm)> Thresh)
subG=induced_subgraph(Net,membership(comm)%in%icomm)
#icomm=which(sizes(comm)>20)
#subG=induced_subgraph(Net,membership(comm)%in%icomm)
#subcomm=wt[icomm]
#Net=subG
#comm=wt
Net=subG

V(Net)$size = 4
nc=length(icomm)
newvals=c(1:nc)
V(Net)$color=newvals[as.factor(V(Net)$comm)]
plot(Net, layout=layout.fruchterman.reingold, vertex.color=V(Net)$color, vertex.size = V(Net)$size, vertex.label = NA, edge.arrow.size = 0.5)

########to store the plotting
#png(filename = "net.png", height = 800, width = 800)
#plot(...)
#dev.off()

