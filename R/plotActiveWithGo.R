source('queryGo.R')
library(igraph)
netdir='../OscopeData/dataAdjacencyMatrixe02/'
#filenet='WaterfallbootstrapSequential_CPUDPSF_TF_2000Boot_data_s45res_C28G8442_edgeNetwork.csv'
filenet='SVZbootstrap_CPUDPSF_TF_2000Boot_data_aNSC_C32G7883_edgeNetwork.csv'
Thresh=8
Dircom='../Communities/Networks02/'
#file='Go_W45fastgreedy.txt'
file='Safastgreedy.txt'
titlestr='SVZ Active Communities (size>=8), FastGreedy algorithm'
ont='BP'
Nlimit=2
plim=0.05

filename=paste(netdir,filenet,sep="")
data=read.csv(file=filename, header= TRUE,row.names=1)
G=graph_from_data_frame(data,directed=FALSE)

#fastgreedy.community
fg=fastgreedy.community(G)

comm=fg
Net=G
V(Net)$comm=membership(comm)
icomm=which(sizes(comm)>= Thresh)
subG=induced_subgraph(Net,membership(comm)%in%icomm)
#icomm=which(sizes(comm)>20)
#subG=induced_subgraph(Net,membership(comm)%in%icomm)
#subcomm=wt[icomm]
#Net=subG
#comm=wt
Net=subG

V(Net)$size = 4
nc=length(icomm)
colors<- rainbow(nc, alpha=.5)
#c(1:nc)####indice valori colori communiti, con questi vanno mapati i nomi
V(Net)$comm=paste('comm',V(Net)$comm,sep='')
ff=as.factor(V(Net)$comm)
V(Net)$color=colors[ff]

plot(Net, layout=layout.fruchterman.reingold, vertex.color=V(Net)$color, vertex.size = V(Net)$size, vertex.label = NA, edge.arrow.size = 0.5)
title(titlestr)

#nn=names(comm[icomm])
#ids=unique(V(Net)$comm)#paste('comm',nn,sep='')

out=queryGo(Dircom,file,ont,Nlimit,plim)
aaa=out[[2]]
labels=aaa[levels(ff),2]

#legend(x=-1.5, y=-0.9, legend=labels, pch=21, pt.bg= colors, pt.cex=2, cex=.8, bty="n", ncol=2,col=colors)
legend(x=-1.5, y=-1, legend=labels, pch=21, pt.bg= colors, pt.cex=2, cex=.8, bty="n", ncol=2,col=colors)






#################
Net=subG

V(Net)$size = 1
E(Net)$count = 1
nc=length(icomm)
newvals=c(1:nc)
V(Net)$comm=newvals[as.factor(V(Net)$comm)]
comm.graph <- contract.vertices(Net, V(Net)$comm, vertex.attr.comb=list(size="sum", "ignore"))
comm.graph <- simplify(comm.graph, remove.loops=TRUE, edge.attr.comb=list(count="sum", "ignore"))
# comm.graph <- contract.vertices(Net, comm$membership, vertex.attr.comb=list(size="sum", "ignore"))
 # comm.graph <- simplify(comm.graph, remove.loops=TRUE, edge.attr.comb=list(count="sum", "ignore"))
 # ###E(comm.graph)$count is the number of edges between communities linked by the edge
 # ###V(comm.graph)$size number of vertex in each community
 M=max(E(comm.graph)$count)
 
 
#karate <- graph.famous("Zachary")
#wckarate <- walktrap.community(karate) #any algorithm

#comm=wckarate
#Net=karate


####play around with scale for Vertex size to improove visualization.

####This is nice for no more then 20 Communities
vectsizes=V(comm.graph)$size
scale=50
#colbar <- rainbow(nc)


