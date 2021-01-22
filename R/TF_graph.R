# Generate TF graph
#' @export
#' @import igraph
TF_graph = function(nTarget=10,nTF=3){
  np=(nTarget+1)*nTF
  adjmat=matrix(0,nrow=np,ncol=np)
  for(j in 1:nTF){
    mj=(nTarget+1)*(j-1)+1
    mj2=(nTarget+1)*j
    adjmat[mj,(mj+1):mj2]=1
    adjmat[(mj+1):mj2,mj]=1
  }
  graph.gen=graph.adjacency(adjmat,mode="undirected")
  return(graph.gen)
}
