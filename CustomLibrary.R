community.significance.test <- function(graph, Nodeids, ...) {
  if (is.directed(graph)) stop("The graph is not undirected")
  subgraph <- induced.subgraph(graph, Nodeids)
  in.degrees <- degree(subgraph)
  out.degrees <- degree(graph, Nodeids) - in.degrees
  res=wilcox.test(in.degrees, out.degrees, alternative = "greater",paired=FALSE,exact=FALSE)
  return(res)
}
