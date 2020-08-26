######create a summary network with nodes given by the communities

####find intra-cluster density:
library(igraph)
karate <- graph.famous("Zachary")
wckarate <- walktrap.community(karate) #any algorithm

comm=wckarate
Net=karate
comm=wt
Net=G
######compute the number of edges within each community
intra=sapply(unique(membership(comm)), function(g) {
    subg1<-induced.subgraph(Net, which(membership(comm)==g)) #membership id differs for each cluster
    ecount(subg1)/ecount(Net)
})
####in order to get the edges between the communities, 
####first we get all combinations of communities

cs <- data.frame(combn(unique(membership(comm)),2))
######compute the edges between each couple of communities
cx <- sapply(cs, function(x) {
    es<-E(Net)[V(Net)[membership(comm)==x[1]] %--% 
              V(Net)[membership(comm)==x[2]]]    
    length(es)
})
between=cbind(t(cs),cx)
between=between[(between[,3]!=0),]
colnames(between)=c('C1','C2','edges')
#create a Graph with communities per nodes and number of edges "between" communities per edges
names=unique(membership(comm))
nc=length(names)
sizesc=NULL
for (i in c(1:nc)){
	sizesc[i]=length(which(membership(comm)==names[i]))
}
communities=data.frame(name=names, intraedge=intra, size=sizesc)
between=as.dataframe(between)
g <- graph_from_data_frame(between, directed=FALSE, vertices=communities)

######to set the color based on a vertex attribute, if na is the number of different values in the attrbute, define a vector of na colors
####if we choose name as example of vertex attribute
colbar <- rainbow(nc)
V(g)$color <- colbar



#colrs <- c("gray50", "tomato", "gold", "blue", "yellow")

#V(g)$color <- colrs[V(g)$name]
 ###plot the edges width according to the number of edges etween communities
 E(g)$width <- E(g)$edges
plot(g,vertex.label.dist=1.5)
#





