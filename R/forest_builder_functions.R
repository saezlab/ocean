#' buildtree
#'
#' @param root  ipsum...
#' @param graph ipsum...
#' @return ipsum...
#' @importFrom igraph neighbors
buildTree <- function(root, graph)
{
  #print(root)
  
  ##initialise variables
  up_down <- list(0)
  k <- 1
  
  #first, upstream (in) and then downstream (out) of the metabolite
  for (mymode in c("in","out"))
  {
    ##initialise
    enzyme_layers <- list(0)
    j <- 1
    
    ##first we get the first neighbors of the metabolite. These are the enzymes that directly metabolise this metabolite. We also intilise all_enzymes, that will help keeping track of enzymes that have already been stored.
    all_enzymes <- neighbors(graph, root, mode = mymode)
    all_enzymes <- unique(names(all_enzymes))
    #print(all_enzymes)
    enzymes <- all_enzymes
    enzyme_layers[[j]] <- enzymes
    
    ##then we repeat a loop that get the next layer of enzymes and keep going as long as new enzymes are found in the next layer.
    while (!is.null(enzymes) & length(enzymes) > 0)
    {
      ##First, the next layer is metabolites. We just use it to reach the next layer of enzymes
      metabs <- list(0)
      for (i in 1:length(enzymes))
      {
        #print(enzymes[i])
        metabs[[i]] <- neighbors(graph, enzymes[i], mode = mymode)
        
      }
      metabs <- unique(names(unlist(metabs)))
      
      #print(metabs)
      if (!is.null(metabs))
      {
        enzymes <- list(0)
        for (i in 1:length(metabs))
        {
          enzymes[[i]] <- neighbors(graph, metabs[i], mode = mymode)
        }
        enzymes <- unique(names(unlist(enzymes)))
        
        #We remove from the current layer the enzyme that already belong to previous layers
        enzymes <- enzymes[!enzymes %in% all_enzymes]
        all_enzymes <- unique(c(all_enzymes,enzymes))
        
        ##go to next step and store the current layer
        j <- j+1
        
        #print(enzymes)
        if (!is.null(enzymes) & length(enzymes) != 0)
        {
          enzyme_layers[[j]] <- enzymes
        }
      }
      else
      {
        enzymes <- NULL
      }
      ##then we get the neighbors of the metabolites. This way we go one layer of metabolic enzymes further.
    }
    ##we store upstream enzymes, then update k and repeat the process for downstream enzymes
    up_down[[k]] <- enzyme_layers
    k <- k+1
  }
  names(up_down) <- c("upstream","downstream")
  return(up_down)
}

#' forestMaker
#'
#' @param molecule_names  ipsum...
#' @param reaction_network ipsum...
#' @return ipsum...
#' @importFrom igraph graph_from_data_frame
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#' @export
forestMaker <- function(molecule_names, reaction_network, branch_length = c(1,1))
{
  graph <- graph_from_data_frame(d = reaction_network[, c(1, 2)], directed = TRUE)
  
  ncores <- min(c(1000,detectCores()-1))
  
  print(ncores)
  
  breaks <- seq(1,length(molecule_names), ncores)
  
  print(breaks)
  
  nruns <- length(breaks)
  
  print(nruns)
  
  splitted_forest <- list()
  if (nruns > 1)
  {
    for (i in 1:(nruns-1))
    {
      print(i)
      molecule_names_run_i <- molecule_names[breaks[i]:(breaks[i+1]-1)]
      print(molecule_names_run_i)
      splitted_forest[[i]] <- mclapply(molecule_names_run_i, buildTree, graph, mc.cores = ncores)
    }
    
    molecule_names_run_i <- molecule_names[breaks[nruns]:length(molecule_names)]
    print(molecule_names_run_i)
    splitted_forest[[nruns]] <- mclapply(molecule_names_run_i, buildTree, graph, mc.cores = ncores)
    
    forest <- do.call(c, splitted_forest)
  }
  else
  {
    print(molecule_names)
    forest <- mclapply(molecule_names, buildTree, graph, mc.cores = ncores)
  }
  
  names(forest) <- molecule_names
  
  forest <- lapply(forest, function(x){
    if(length(x[[1]]) > branch_length[1] | length(x[[2]]) > branch_length[2])
    {
      return(x)
    } else
    {
      return(NA)
    }
  })
  
  forest <- forest[!is.na(forest)]
  
  return(forest)
}