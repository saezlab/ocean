#'\code{model_to_pathway_sif}
#'
#' This function generate a subreaction network in sif format from the recon_redhuman model for a given set of pathways
#'
#' @param pathway_to_keep  ipsum...
#' @param pathways ipsum...
#' @param rxns ipsum...
#' @param rules ipsum...
#' @param minmax ipsum...
#' @param stochio ipsum...
#' @param metab ipsum...
#' @param mapping ipsum...
#' @param recon1_metabolites ipsum...
#' @return ipsum...
#' @export

model_to_pathway_sif <- function(pathway_to_keep, 
                                 pathways = recon2_redhuman$pathway,  #redHuman_models folder
                                 rxns = recon2_redhuman$rxns,
                                 rules = recon2_redhuman$rules,
                                 minmax = recon2_redhuman$minmax,
                                 stochio = recon2_redhuman$stochio,
                                 metab = recon2_redhuman$metab,
                                 mapping = recon2_redhuman$gene_mapping,
                                 recon1_metabolites = recon2_redhuman$metab_mapping){
  ####
  
  index_to_keep <- sapply(pathway_to_keep, function(x, pathways){
    which(pathways == x)
  }, pathways = pathways, simplify = T)
  
  index_to_keep <- unique(unlist(index_to_keep))
  
  ####
  
  rxns <- rxns[index_to_keep,, drop = F]
  
  rules <- rules[index_to_keep,, drop = F]
  
  minmax <- minmax[index_to_keep,, drop = F]
  
  pathways <- pathways[index_to_keep,, drop = F]
  
  stochio <- stochio[,index_to_keep,drop = F]
  
  ####
  
  print("Building gprs...")
  gprs <- build_grps(rxns, rules, minmax, pathways)
  ####
  
  print("Building reaction list...")
  reactions_list <- build_reactions_list(gprs, stochio, metab, minmax)
  ###
  
  print("Building reaction data frame...")
  reactions_df <- build_reactions_df(reactions_list, mapping)
  ###
  
  # print("Converting network gene IDs...")
  # reactions_df <- translate_network_gene_ids(reactions_df, mapping)

  
  reactions_df <- reactions_df[!grepl("HC0",reactions_df$V1) & !grepl("HC0",reactions_df$V2),]
  reactions_df <- reactions_df[complete.cases(reactions_df),]
  
  ###
  
  print("Converting network metabolite IDs...")
  reactions_df <- translate_network_metab_ids(reactions_df, recon1_metabolites)
  
  names(reactions_df) <- c("source","target")
  
  reactions_df <- unique(reactions_df)
  
  ###
  
  row.names(reactions_df) <- c(1:length(reactions_df[,1]))
  
  if("Urea cycle" %in% pathway_to_keep)
  {
    print("Ataching Urea cycle additional reactions")
    reactions_df <- as.data.frame(rbind(reactions_df,network_supplements[["Urea cycle"]]))
  }
  if("Purine synthesis" %in% pathway_to_keep)
  {
    print("Ataching Purine synthesis additional reactions")
    reactions_df <- as.data.frame(rbind(reactions_df,network_supplements[["Purine synthesis"]]))
  }
  if("Valine, leucine, and isoleucine metabolism" %in% pathway_to_keep)
  {
    print("Ataching BCAA additional reactions")
    reactions_df <- as.data.frame(rbind(reactions_df,network_supplements[["Carnitine synthesis"]])) #add ketoacid to carnitine
    
    reactions_df <- as.data.frame(rbind(reactions_df,network_supplements[["Valine leucine and isoleucine metabolism"]])) #add bcat1 reverse
    
    reactions_df <- reactions_df[reactions_df[,1] != "1629_594_593_1738" & #remove the wrong BCKDH complex and correct it
                                   reactions_df[,2] != "1629_594_593_1738" &
                                   reactions_df[,1] != "1629_1738_594_593" &
                                   reactions_df[,2] != "1629_1738_594_593",]
    
    reactions_df <- reactions_df[!grepl("CRAT.*_reverse",reactions_df[,1]) & #remove reverse crat
                                   !grepl("CRAT.*_reverse",reactions_df[,2]),]
    
    reactions_df <- as.data.frame(rbind(reactions_df,network_supplements[["BCKDH"]])) #add correct BCKDH
  }
  ###
  return(reactions_df)
}

#'\code{remove_cofactors}
#'
#' This function removes cofactors from reaction networks generate by the model_to_pathway_sif function
#'
#' @param reaction_network  ipsum...
#' @param compound_list KEGGREST coumpound list
#' @return ipsum...
#' @export
remove_cofactors <- function(reaction_network, compound_list = kegg_compounds)
{
  compounds <- unique(c(reaction_network$source, reaction_network$target))
  
  compounds <- compounds[grepl("cpd[:]",compounds)]
  
  compounds <- as.data.frame(compounds)
  names(compounds)[1] <- "compounds_compartment"
  compounds$compound <- gsub("_.*","",compounds$compounds_compartment)
  
  compound_vec <- unique(compounds$compound)
  
  bad_kegg_compounds <- list()
  for (compound in compound_list)
  {
    
    brite <- compound$BRITE
    compound$NAME <- gsub(";","",compound$NAME)
    if(!is.null(compound$ATOM))
    {
      brite <- gsub("[ ]+","",brite)
      if(("Cofactors" %in% brite | 
         "Nucleotides" %in% brite | 
         "CO2" %in% compound$NAME | 
         as.numeric(compound$ATOM[1]) <= 3 | 
         "ITP" %in% compound$NAME | 
         "IDP" %in% compound$NAME | 
         "NADH" %in% compound$NAME | 
         "NADPH" %in% compound$NAME | 
         "FADH2" %in% compound$NAME) &
         !"SAM" %in% compound$NAME)
      {
        bad_kegg_compounds[[compound$ENTRY[1]]] <- compound
      }
    }
  }
  
  bad_entries <- c()
  for(compound in bad_kegg_compounds)
  {
    bad_entries <- c(bad_entries,compound$ENTRY[1])
  }
  
  bad_entries <- paste("cpd:",bad_entries,sep = "")
  
  bad_compounds <- compounds[compounds$compound %in% bad_entries,]
  
  ##manually add diphosphate to be removed from network
  bad_compounds <- as.data.frame(rbind(bad_compounds,c("cpd:C00013_c","cpd:C00013")))
  bad_compounds <- as.data.frame(rbind(bad_compounds,c("cpd:C00013_m","cpd:C00013")))
  bad_compounds <- as.data.frame(rbind(bad_compounds,c("cpd:C00009_c","cpd:C00009")))
  bad_compounds <- as.data.frame(rbind(bad_compounds,c("cpd:C00009_m","cpd:C00009")))
  bad_compounds <- as.data.frame(rbind(bad_compounds,c("cpd:C00009_r","cpd:C00009")))
  
  reaction_network_no_cofact <- reaction_network[
    !(reaction_network$source %in% bad_compounds$compounds_compartment) & 
      !(reaction_network$target %in% bad_compounds$compounds_compartment),
  ]
  
  reaction_network_no_cofact <- 
    reaction_network_no_cofact[
      complete.cases(reaction_network_no_cofact),
    ]
  
  nodes <- unique(c(reaction_network_no_cofact$source ,
                    reaction_network_no_cofact$target ))
  
  node_attributes <- as.data.frame(matrix(NA,length(nodes),2))
  node_attributes[,1] <- nodes
  node_attributes[,2] <- "reaction_enzyme"
  node_attributes[node_attributes[,1] %in% compounds$compounds_compartment,2] <- "metabolite"
  
  return(list("reaction_network" = reaction_network_no_cofact, "attributes" = node_attributes))
}

#'\code{compress_transporters}
#'
#' This function removes cofactors from reaction networks generate by the model_to_pathway_sif function
#'
#' @param sub_network_nocofact  ipsum...
#' @return ipsum...
#' @export
compress_transporters <- function (sub_network_nocofact) 
{
  test_1 <- paste(gsub("_[cxrnme]$", "", sub_network_nocofact$reaction_network$source), 
                  gsub("_[cxrnme]$", "", sub_network_nocofact$reaction_network$target), 
                  sep = "_")
  test_2 <- paste(gsub("_[cxrnme]$", "", sub_network_nocofact$reaction_network$target), 
                  gsub("_[cxrnme]$", "", sub_network_nocofact$reaction_network$source), 
                  sep = "_")
  test_1[test_1 %in% test_2]
  transporters <- unique(c(sub_network_nocofact$reaction_network[test_1 %in% 
                                                                   test_2, 1], sub_network_nocofact$reaction_network[test_1 %in% 
                                                                                                                       test_2, 2]))
  transporters <- transporters[!grepl("_[cxrnme]$", transporters)]
  for (i in 1:2) {
    sub_network_nocofact$reaction_network[, i] <- sapply(sub_network_nocofact$reaction_network[, 
                                                                                               i], function(x, transporters) {
                                                                                                 if (x %in% transporters) {
                                                                                                   return("transporter")
                                                                                                 }
                                                                                                 else {
                                                                                                   return(x)
                                                                                                 }
                                                                                               }, transporters = transporters)
  }
  sub_network_nocofact$reaction_network$edgeId <- paste(sub_network_nocofact$reaction_network$source, 
                                                        sub_network_nocofact$reaction_network$target, sep = "_")
  sub_network_nocofact$reaction_network <- sub_network_nocofact$reaction_network[!duplicated(sub_network_nocofact$reaction_network$edgeId), 
  ]
  edge_transporter <- sub_network_nocofact$reaction_network$edgeId
  edge_transporter <- edge_transporter[grepl("transporter", 
                                             edge_transporter)]
  groups <- rep(0, length(edge_transporter))
  metab_memory <- c()
  for (i in 1:length(edge_transporter)) {
    if (grepl("^cpd:", edge_transporter[i]) | grepl("_[cxrnme]_", 
                                                    edge_transporter[i])) {
      metab <- gsub("_.*", "", edge_transporter[i])
      if(!metab %in% metab_memory)
      {
        metab_memory[i] <- metab
        for (j in 1:length(edge_transporter)) {
          if (grepl(metab, edge_transporter[j]))
          {
            groups[j] <- i
          }
        }
      }
    }
  }
  names(groups) <- edge_transporter
  network <- sub_network_nocofact$reaction_network
  for (i in 1:2) {
    for (j in 1:length(network[, i])) {
      if (network[j, i] == "transporter") {
        print(network[j, 3])
        network[j, i] <- paste(network[j, i], groups[network[j, 
                                                             3]], sep = ">")
      }
    }
  }
  network <- network[!grepl("^SLC", network$edgeId) & !grepl("_SLC", 
                                                             network$edgeId), ]
  sub_network_nocofact$reaction_network <- network[, -3]
  sub_network_nocofact$attributes <- sub_network_nocofact$attributes[sub_network_nocofact$attributes[, 
                                                                                                     1] %in% network[, 1] | sub_network_nocofact$attributes[, 
                                                                                                                                                            1] %in% network[, 2], ]
  new_transporters <- unique(c(network[, 1], network[, 2]))
  new_transporters <- new_transporters[grepl("transporter", 
                                             new_transporters)]
  new_transporters <- as.data.frame(cbind(new_transporters, 
                                          rep("transporter", length(new_transporters))))
  names(new_transporters) <- c("V1", "V2")
  sub_network_nocofact$attributes <- as.data.frame(rbind(sub_network_nocofact$attributes, 
                                                         new_transporters))
  return(sub_network_nocofact)
}

#'\code{split_transaminases}
#'
#' This function removes cofactors from reaction networks generate by the model_to_pathway_sif function
#'
#' @param sub_network_nocofact  ipsum...
#' @return ipsum...
#' @export
split_transaminases <- function(sub_network_nocofact)
{
  sub_net <- sub_network_nocofact$reaction_network
  
  transaminases_net <- sub_net[grepl("C0002[56]",sub_net$source) |
                                 grepl("C0002[56]",sub_net$target) |
                                 grepl("C00064",sub_net$source) |
                                 grepl("C00064",sub_net$target),]
  
  
  sub_net_no_transaminase <- sub_net[!(grepl("C0002[56]",sub_net$source) |
                                         grepl("C0002[56]",sub_net$target) |
                                         grepl("C00064",sub_net$source) |
                                         grepl("C00064",sub_net$target)),]
  
  transaminases_potential <- unique(c(transaminases_net$source,transaminases_net$target))
  transaminases_potential <- transaminases_potential[!grepl("cpd:",transaminases_potential)]
  
  transaminases_net <- sub_net[sub_net$source %in% transaminases_potential |
                                 sub_net$target %in% transaminases_potential,]
  
  splitted_transaminase <- sapply(transaminases_potential, function(transaminase, transaminases_net){
    sub_net_transaminase <- transaminases_net[transaminase == transaminases_net$source |
                                                transaminase == transaminases_net$target,]
    
    columns <- c(1,2)
    
    elements <- unique(c(sub_net_transaminase[,1],sub_net_transaminase[,2]))
    
    if(sum(grepl("cpd:C00025_.",elements)) & sum(grepl("cpd:C00026_.",elements)))
    {
      for(i in 1:length(sub_net_transaminase[,1]))
      {
        for(j in 1:2)
        {
          if(grepl("C0002[56]",sub_net_transaminase[i,j]))
          {
            sub_net_transaminase[i,columns[-j]] <- paste(sub_net_transaminase[i,columns[-j]],"_gluakg",sep = "")
          }
        }
      }
    }
    
    if(sum(grepl("cpd:C00025_.",elements)) & sum(grepl("cpd:C00064_.",elements)))
    {
      for(i in 1:length(sub_net_transaminase[,1]))
      {
        for(j in 1:2)
        {
          if(grepl("C00025",sub_net_transaminase[i,j]) | grepl("C00064",sub_net_transaminase[i,j]))
          {
            sub_net_transaminase[i,columns[-j]] <- paste(sub_net_transaminase[i,columns[-j]],"_glugln",sep = "")
          }
        }
      }
    }
    
    return(sub_net_transaminase)
  }, transaminases_net = transaminases_net, USE.NAMES = T, simplify = F)
  
  splitted_transaminase <- as.data.frame(do.call(rbind,splitted_transaminase))
  
  sub_network_nocofact$reaction_network <- as.data.frame(rbind(sub_net_no_transaminase, splitted_transaminase))
  
  sub_network_nocofact$reaction_network <- unique(sub_network_nocofact$reaction_network)
  
  network_attributes <- as.data.frame(cbind(unique(c(sub_network_nocofact$reaction_network[,1],sub_network_nocofact$reaction_network[,2])),NA))
  network_attributes[,2] <- ifelse(grepl("cpd:",network_attributes[,1]), "metabolite","reaction_enzyme")
  names(network_attributes) <- names(sub_network_nocofact$attributes)
  
  sub_network_nocofact$attributes <- network_attributes
  
  return(sub_network_nocofact)
}

#'\code{nitrogen_tracking}
#'
#' This function cross link transamination metabolites to track amine groups
#'
#' @param sub_network_nocofact  ipsum...
#' @return ipsum...
#' @export
nitrogen_tracking <- function(sub_network_nocofact)
{
  sub_net <- sub_network_nocofact$reaction_network
  
  transaminases_net <- sub_net[grepl("C0002[56]",sub_net$source) |
                                 grepl("C0002[56]",sub_net$target),]
  
  
  
  
  transaminases_potential <- unique(c(transaminases_net$source,transaminases_net$target))
  transaminases_potential <- transaminases_potential[!grepl("cpd:",transaminases_potential)]
  
  transaminases_net <- sub_net[sub_net$source %in% transaminases_potential |
                                 sub_net$target %in% transaminases_potential,]
  
  sub_net_no_transaminase <- sub_net[!(sub_net$source %in% transaminases_potential |
                                         sub_net$target %in% transaminases_potential),]
  
  
  splitted_transaminase <- sapply(transaminases_potential, function(transaminase, transaminases_net){
    sub_net_transaminase <- transaminases_net[transaminase == transaminases_net$source |
                                                transaminase == transaminases_net$target,]
    
    columns <- c(1,2)
    
    elements <- unique(c(sub_net_transaminase[,1],sub_net_transaminase[,2]))
    
    if(sum(grepl("cpd:C00025_.",elements)) & sum(grepl("cpd:C00026_.",elements)))
    {
      if(sum(grepl("C00025",sub_net_transaminase[,1])) > 0)
      {
        for(i in 1:length(sub_net_transaminase[,1]))
        {
          if(grepl("C00025",sub_net_transaminase[i,1]))
          {
            sub_net_transaminase[i,2] <- paste(sub_net_transaminase[i,2],"_nitro",sep = "")
          } 
          if(!grepl("C00026",sub_net_transaminase[i,2]) & grepl("cpd:C",sub_net_transaminase[i,2]))
          {
            sub_net_transaminase[i,1] <- paste(sub_net_transaminase[i,1],"_nitro",sep = "")
          } 
        }
      }else
      {
        for(i in 1:length(sub_net_transaminase[,1]))
        {
          if(grepl("C00025",sub_net_transaminase[i,2]))
          {
            sub_net_transaminase[i,1] <- paste(sub_net_transaminase[i,1],"_nitro",sep = "")
          } 
          if(!grepl("C00026",sub_net_transaminase[i,1]) & grepl("cpd:C",sub_net_transaminase[i,1]))
          {
            sub_net_transaminase[i,2] <- paste(sub_net_transaminase[i,2],"_nitro",sep = "")
          } 
        }
      }
    }
    return(sub_net_transaminase)
  }, transaminases_net = transaminases_net, USE.NAMES = T, simplify = F)
  
  splitted_transaminase <- as.data.frame(do.call(rbind,splitted_transaminase))
  
  sub_network_nocofact$reaction_network <- as.data.frame(rbind(sub_net_no_transaminase, splitted_transaminase))
  
  sub_network_nocofact$reaction_network <- unique(sub_network_nocofact$reaction_network)
  
  network_attributes <- as.data.frame(cbind(unique(c(sub_network_nocofact$reaction_network[,1],sub_network_nocofact$reaction_network[,2])),NA))
  network_attributes[,2] <- ifelse(grepl("cpd:",network_attributes[,1]), "metabolite","reaction_enzyme")
  names(network_attributes) <- names(sub_network_nocofact$attributes)
  
  sub_network_nocofact$attributes <- network_attributes
  
  return(sub_network_nocofact)
}
