
library(igraph)
library(gridExtra)
library(VennDiagram)

######Venn Diagrams
#https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html

####GRAPH filename Very Stringent case

filename1='../OscopeData/dataAdjacencyMatrixe02/SVZbootstrap_CPUDPSF_TF_2000Boot_data_qNSC_C60G7883_edgeNetwork.csv'
filename2='../OscopeData/dataAdjacencyMatrixe02/SVZbootstrap_CPUDPSF_TF_2000Boot_data_aNSC_C32G7883_edgeNetwork.csv'

data1=read.csv(file=filename1, header= TRUE,row.names=1)
data2=read.csv(file=filename2, header= TRUE,row.names=1)

G1=graph_from_data_frame(data1,directed=FALSE)
G2=graph_from_data_frame(data2,directed=FALSE)
g1name=get.vertex.attribute(G1)$name
g2name=get.vertex.attribute(G2)$name
g1g2=intersect(g1name, g2name)
ncomm=length(g1g2)
print(ncomm)
print(length(g1name))
print(length(g2name))

#grid.newpage()
venn.plot.df1=draw.pairwise.venn(area1 = length(g1name), area2 = length(g2name), cross.area = ncomm, category = c("(a) quiscent", 
    "(a) active"), lty = rep("blank",2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))



####GRAPH filename Almost Stringent case

filename1='../OscopeData/dataAdjacencyMatrixePiDiv9/SVZbootstrap_CPUDPSF_TF_2000Boot_data_qNSC_C60G7883_edgeNetwork.csv'
filename2='../OscopeData/dataAdjacencyMatrixePiDiv9/SVZbootstrap_CPUDPSF_TF_2000Boot_data_aNSC_C32G7883_edgeNetwork.csv'

data1=read.csv(file=filename1, header= TRUE,row.names=1)
data2=read.csv(file=filename2, header= TRUE,row.names=1)

G1=graph_from_data_frame(data1,directed=FALSE)
G2=graph_from_data_frame(data2,directed=FALSE)
g1name=get.vertex.attribute(G1)$name
g2name=get.vertex.attribute(G2)$name
g1g2=intersect(g1name, g2name)
ncomm=length(g1g2)
print(ncomm)
print(length(g1name))
print(length(g2name))

#grid.newpage()
venn.plot.df2=draw.pairwise.venn(area1 = length(g1name), area2 = length(g2name), cross.area = ncomm, category = c("(b) quiscent", 
    "(b) active"), lty = rep("blank",2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))
filename1='../OscopeData/FirstData/SVZbootstrap_CPUDPSF_TF_2000Boot_data_qNSC_C60G7883_edgeNetwork.csv'
filename2='../OscopeData/FirstData/SVZbootstrap_CPUDPSF_TF_2000Boot_data_aNSC_C32G7883_edgeNetwork.csv'


####GRAPH filename non stringent case

data1=read.csv(file=filename1, header= TRUE,row.names=1)
data2=read.csv(file=filename2, header= TRUE,row.names=1)

G1=graph_from_data_frame(data1,directed=FALSE)
G2=graph_from_data_frame(data2,directed=FALSE)
g1name=get.vertex.attribute(G1)$name
g2name=get.vertex.attribute(G2)$name
g1g2=intersect(g1name, g2name)
ncomm=length(g1g2)
print(ncomm)
print(length(g1name))
print(length(g2name))

grid.newpage()
venn.plot.df3=draw.pairwise.venn(area1 = length(g1name), area2 = length(g2name), cross.area = ncomm, category = c("(c) quiscent", 
    "(c) active"), lty = rep("blank",2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(10,20), cat.dist = rep(0.025, 2))



venn.plot=grid.arrange(gTree(children=venn.plot.df1,name='primo'), gTree(children=venn.plot.df2),
             gTree(children=venn.plot.df3), ncol=3,top="Quiescent-Active Cooscillating Genes", bottom="Three Scenario: (a) Very Stringent (c) Partially Stringent (b) Non Stringent")
             
#tiff(filename = "Global_Venn_diagram.tiff", compression = "lzw");
tiff(filename = "Global_Venn_diagramSVZ.tiff",width=780,height=480,pointsize=20);

#jpeg
grid.draw(venn.plot);
dev.off();


