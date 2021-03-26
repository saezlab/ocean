

#' make_discrete_palette
#'
#' @param numeric_vector  a metabolic-enzyme set collections, 
#' as returned by the condense_metabolite_set function
#' @param rbrewer_plalette_name numeric, the index of the reaction_set_list_merged corresponding to
#' the desired penalty
#' @return ipsum lorem...
#' @importFrom scales rescale
#' @importFrom RColorBrewer brewer.pal
make_discrete_palette <- function(numeric_vector, rbrewer_plalette_name)
{
  mirrored <- c(numeric_vector, -numeric_vector)
  
  mirrored <- scales::rescale(mirrored, to = c(1,11))
  
  binned <- as.numeric(cut(mirrored, breaks = 1:11, include.lowest = T))
  
  my_palette <- RColorBrewer::brewer.pal(n = max(binned), name = rbrewer_plalette_name)
  
  colors <- sapply(binned, function(x, my_palette){
    my_palette[x]
  }, my_palette = my_palette)
  
  colors <- colors[1:length(numeric_vector)]
  
  return(colors)
}

#' plot_reaction_network
#'
#' 
#'
#' @param network_and_attributes  lorem ipsum...
#' @param t_table lorem ipsum...
#' @param scores_df lorem ipsum...
#' @param column_index lorem ipsum...
#' @param rbrewer_plalette_name lorem ipsum...
#' @return lorem ipsum...
#' @export
#' @importFrom visNetwork visNetwork
#' @importFrom visNetwork visOptions
plot_reaction_network <- function(network_and_attributes, t_table, scores_df, column_index, rbrewer_plalette_name = "RdYlGn", vis.height = 700, vis.degree = 2)
{
  
  edges <- network_and_attributes[[1]]
  
  edges <- edges[edges[,1] %in% scores_df[,1] | edges[,2] %in% scores_df[,1],]
  
  nodes <- network_and_attributes[[2]]
  
  nodes <- nodes[nodes[,1] %in% edges[,1] | nodes[,1] %in% edges[,2],]
  
  names(nodes) <- c("id","molecule_type")
  
  metab_stat <- t_table[,c(1,column_index+1)]
  names(metab_stat) <- c("id","metab_stat")
  
  enzyme_score <- scores_df[,c(1,column_index+1)]
  names(enzyme_score) <- c("id","enzyme_score")
  
  nodes <- merge(nodes, metab_stat, all.x = T)
  
  nodes <- merge(nodes, enzyme_score, all.x = T)
  
  nodes$metab_stat[is.na(nodes$metab_stat)] <- 0
  
  nodes$enzyme_score[is.na(nodes$enzyme_score)] <- 0
  
  nodes$stat <- nodes$metab_stat + nodes$enzyme_score
  
  nodes$color <- make_discrete_palette(nodes$stat, rbrewer_plalette_name)
  
  nodes$color <- ifelse(nodes$stat == 0, "#C1C1C1", nodes$color)
  
  nodes$label <- nodes$id
  
  nodes$shape <- ifelse(nodes$molecule_type == "reaction_enzyme", "square", "triangle")
  nodes$shape <- ifelse(grepl("transporter",nodes$id), "circle", nodes$shape)
  
  nodes$label <- gsub("^transporter.*","",nodes$label)
  
  names(edges) <- c("from","to")
  
  edges$arrows <- "to"
  
  mapping_vec <- mapping_table$metab
  names(mapping_vec) <- mapping_table$KEGG
  
  nodes$label <- sapply(nodes$label, function(x,mapping_vec){
    
    if(grepl("cpd:",x))
    {
      suffixe <- stringr::str_extract(x, "_.$")
      x <- gsub("_.$","",x)
      x <- gsub("cpd:","",x)
      if(x %in% names(mapping_vec))
      {
        x <- mapping_vec[x]
      }
      
      if(!is.na(suffixe))
      {
        x <- paste0(x, suffixe)
      }
      
      return(x)
    } else
    {
      return(x)
    }
  }, mapping_vec = mapping_vec, simplify = F, USE.NAMES = F)
  
  visNetwork::visNetwork(nodes = nodes, edges = edges, width = "100%", height = vis.height) %>% 
    visOptions(nodesIdSelection = TRUE,
               highlightNearest = list(enabled = T, degree = vis.degree, hover = T))
}