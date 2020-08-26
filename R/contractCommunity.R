# Use the contract.vertices function with the membership vector that the community detection method provides, followed by simplify. In particular:
	# 1.	Assign a numeric vertex attribute with a value of 1 to each vertex as follows: V(g)$size = 1 
	# 2.	Assign a numeric edge attribute with a value of 1 to each edge as follows: E(g)$count = 1 
	# 3.	Contract the communities into vertices as follows: comm.graph <- contract.vertices(g, wc$membership, vertex.attr.comb=list(size="sum", "ignore")); basically this specifies that the size attribute of the vertices being contracted should be summed and every other vertex attribute should be ignored. (See ?attribute.combination in R for more details). This call contracts the vertices but leaves the original edges so you now have as many edges between the vertices as there were in the original graph between the communities. 
	# 4.	Collapse the multiple edges as follows: comm.graph <- simplify(comm.graph, remove.loops=FALSE, edge.attr.comb=list(count="sum", "ignore")). 
# You now have a graph named comm.graph where the vertices represent the communities of the original graph, the size vertex attribute corresponds to the number of vertices in each community in the original graph, and the count edge attribute corresponds to the number of edges between communities in the original graph.

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

######add here the name of the GO group they are mainly in :-)
V(comm.graph)$name <- letters[1:vcount(comm.graph)]
###############Set the best colors range you prefer
V(comm.graph)$color <- c(1:vcount(comm.graph))
V(comm.graph)$size <- vectsizes/max(vectsizes)*scale
E(comm.graph)$weight <-E(comm.graph)$count
E(comm.graph)$length <-M-E(comm.graph)$count#1/E(comm.graph)$count
plot(comm.graph,layout=layout.fruchterman.reingold,vertex.label.dist=1)
#plot(comm.graph,layout= layout_as_star)
#tkplot(comm.graph)


#colrs <- c("gray50", "tomato", "gold", "blue", "yellow")

#V(g)$color <- colrs[V(g)$name]
 ###plot the edges width according to the number of edges etween communities
 #
