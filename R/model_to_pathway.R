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
                                 pathways = recon2_redhuman$pathway, 
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
  
  stochio <- stochio[,index_to_keep]
  
  ####
  
  gprs <- as.data.frame(cbind(rxns,rules))
  
  names(gprs) <- c("Reactions","Genes")
  
  gprs$Genes <- gsub("[()]","", gprs$Genes)
  
  gprs$Genes <- gsub(" or ",";",gprs$Genes)
  gprs$Genes <- gsub(" and ","_",gprs$Genes)
  gprs$Genes <- gsub("[.][0-9]","",gprs$Genes)
  
  gprs$direction <- minmax$direction
  
  gprs$pathway <- pathways
  
  ####
  
  recon1_reactions <- list()
  for(i in 1:length(gprs[,1]))
  {
    recon1_reactions[[i]] <- list()
    reactants <- as.numeric(row.names(stochio[stochio[,i] <= -1,]))
    
    if(!identical(reactants, numeric(0)))
    {
      reactants <- metab[reactants,]
      recon1_reactions[[i]]$reactants <- reactants
    }
    else
    {
      recon1_reactions[[i]]$reactants <- NA
    }
    
    
    products <- as.numeric(row.names(stochio[stochio[,i] >= 1,]))
    
    if(!identical(products, numeric(0)))
    {
      products <- metab[products,]
      recon1_reactions[[i]]$products<- products
    }
    else
    {
      recon1_reactions[[i]]$products <- NA
    }
    recon1_reactions[[i]]$direction <- minmax[i,3]
  }
  names(recon1_reactions) <- gprs$Reactions
  
  for(i in 1:length(recon1_reactions))
  {
    reaction <- names(recon1_reactions)[i]
    
    recon1_reactions[[i]]$genes <- gprs[gprs$Reactions == reaction,2]
  }
  
  ###
  
  SLC_vec <- unique(mapping[grepl("SLC[0-9]",mapping$name),1])
  SLC_vec <- c(SLC_vec,
               "340024_57393",
               "59272",
               "9057_6520",
               "6520_56301",
               "56301_6520",
               "9057_6520",
               "8140_6520",
               "6520_23428",
               "6520_8140",
               "9056_6519",
               "6579",
               "23428_6520",
               "11136_6519")
  
  ###
  
  
  reactions_df <- as.data.frame(matrix(NA,0,2)) 
  
  genes_memory <- c()
  z <- 1
  m <- 1
  for(reaction in recon1_reactions)
  {
    if(!is.null(reaction$reactants) &!is.null(reaction$products))
    {
      if(!is.null(reaction) & !is.na(reaction$genes))
      {
        genes <- strsplit(reaction$genes,";")
        genes <- genes[[1]]
      }
      else
      {
        # genes <- reaction$id
        genes <- names(recon1_reactions)[m]
      }
      
      if(reaction$direction == 1)
      {
        reaction$reversible <- FALSE
        reactants <- reaction$reactants
        products <- reaction$products
      }
      else
      {
        if(reaction$direction != 0)
        {
          reactants <- reaction$products
          products <- reaction$reactants
          reaction$reversible <- FALSE
        }
        else
        {
          reaction$reversible <- TRUE
          reactants <- reaction$reactants
          products <- reaction$products
        }
      }
      
      
      reaction_df <- as.data.frame(matrix(NA,0,2))
      
      # print(genes)
      
      for(gene in genes)
      {
        # "6570" "6571" 
        # if(gene == "6523")
        # {
        #   print(reaction)
        # }
        
        if(!(gene %in% genes_memory))
        {
          genes_memory <- c(genes_memory, gene)
        }
        else
        {
          gene <- paste(gene,z,sep = ">")
          genes_memory <- c(genes_memory, gene)
          z <- z+1
        }
        
        if(gsub("[>].*","",gene) %in% SLC_vec)
        {
          new_reaction_df <- as.data.frame(matrix(NA,length(reactants)+length(products), 2))
          gene_i <- c()
          for(i in 1:length(reactants))
          {
            if(length(reactants) == length(products))
            {
              gene_i <- c(gene_i,paste(gene,i,sep = "<"))
            }
            else
            {
              gene_i <- gene
            }
          }
          new_reaction_df[1:length(reactants),2] <- gene_i
          new_reaction_df[(length(reactants)+1):((length(reactants)+1)+length(products)-1),1] <- gene_i
          
          new_reaction_df[1:length(reactants),1] <- reactants
          new_reaction_df[(length(reactants)+1):((length(reactants)+1)+length(products)-1),2] <- products
          
          if(length(reactants) == length(products)) #this should fix the rare cases were the order or products was inverted
          {
            if(sum(gsub("_[cxrnme]","",reactants) == gsub("_[cxrnme]","",products)) != length(reactants))
            {
              new_reaction_df[(length(reactants)+1):((length(reactants)+1)+length(products)-1),2] <- products[length(products):1]
            }
          }
          
          reverse_reaction_df <- new_reaction_df
          reverse_reaction_df[1:length(reactants),2] <- paste(reverse_reaction_df[1:length(reactants),2],"_reverse",sep = "")
          reverse_reaction_df[(length(reactants)+1):((length(reactants)+1)+length(products)-1),1] <- paste(reverse_reaction_df[(length(reactants)+1):((length(reactants)+1)+length(products)-1),1],"_reverse",sep = "")
          names(reverse_reaction_df) <- names(reverse_reaction_df)[c(2,1)]
          
          new_reaction_df <- as.data.frame(rbind(new_reaction_df,reverse_reaction_df))
        }
        else
        {
          new_reaction_df <- as.data.frame(matrix(NA,length(reactants)+length(products), 2))
          new_reaction_df[1:length(reactants),2] <- gene
          new_reaction_df[(length(reactants)+1):((length(reactants)+1)+length(products)-1),1] <- gene
          
          new_reaction_df[1:length(reactants),1] <- reactants
          new_reaction_df[(length(reactants)+1):((length(reactants)+1)+length(products)-1),2] <- products
          
          if(reaction$reversible)
          {
            reverse_reaction_df <- new_reaction_df
            reverse_reaction_df[1:length(reactants),2] <- paste(reverse_reaction_df[1:length(reactants),2],"_reverse",sep = "")
            reverse_reaction_df[(length(reactants)+1):((length(reactants)+1)+length(products)-1),1] <- paste(reverse_reaction_df[(length(reactants)+1):((length(reactants)+1)+length(products)-1),1],"_reverse",sep = "")
            names(reverse_reaction_df) <- names(reverse_reaction_df)[c(2,1)]
            
            new_reaction_df <- as.data.frame(rbind(new_reaction_df,reverse_reaction_df))
          }
        }
        
        reaction_df <- as.data.frame(rbind(reaction_df,new_reaction_df))
      }
      
      reactions_df <- as.data.frame(rbind(reactions_df, reaction_df))
      #print(length(reactions_df[,1]))
    }
    m <- m+1
  }
  
  ###
  
  mapping_vec <- mapping$name
  names(mapping_vec) <- mapping$X1
  
  for (i in 1:2)
  {
    for (j in 1:length(reactions_df[,1]))
    {
      if(grepl(">",reactions_df[j,i]))
      {
        affixe <- gsub(".+[>]","",reactions_df[j,i])
        reactions_df[j,i] <- gsub("[>].*","",reactions_df[j,i])
        has_affixe <- TRUE
      }
      else
      {
        if(grepl("<",reactions_df[j,i]))
        {
          affixe_lower <- gsub(".+[<]","",reactions_df[j,i])
          reactions_df[j,i] <- gsub("[<].*","",reactions_df[j,i])
          has_affixe_lower <- TRUE
        }
        else
        {
          has_affixe_lower <- FALSE
        }
        has_affixe <- FALSE
      }
      if(reactions_df[j,i] %in% names(mapping_vec))
      {
        if(has_affixe)
        {
          reactions_df[j,i] <- paste(mapping_vec[reactions_df[j,i]],affixe,sep = ">")
        }
        else
        {
          if(has_affixe_lower)
          {
            reactions_df[j,i] <- paste(mapping_vec[reactions_df[j,i]],affixe_lower,sep = "<")
          }
          else
          {
            reactions_df[j,i] <- mapping_vec[reactions_df[j,i]]
          }
        }
      }
      if(gsub("_reverse","",reactions_df[j,i]) %in% names(mapping_vec))
      {
        reactions_df[j,i] <- paste(mapping_vec[gsub("_reverse","",reactions_df[j,i])],"_reverse",sep = "")
      }
    }
  }
  
  reactions_df <- reactions_df[!grepl("HC0",reactions_df$V1) & !grepl("HC0",reactions_df$V2),]
  reactions_df <- reactions_df[complete.cases(reactions_df),]
  
  ###
  
  recon1_metabolites_named_vec <- recon1_metabolites[,2]
  names(recon1_metabolites_named_vec) <- recon1_metabolites[,1]
  
  names(recon1_metabolites_named_vec) <- gsub("__","_",names(recon1_metabolites_named_vec))
  
  for (i in 1:length(reactions_df[,1]))
  {
    if (reactions_df[i,1] %in% names(recon1_metabolites_named_vec))
    {
      reactions_df[i,1] <- recon1_metabolites_named_vec[reactions_df[i,1]]
    }
    if (reactions_df[i,2] %in% names(recon1_metabolites_named_vec))
    {
      reactions_df[i,2] <- recon1_metabolites_named_vec[reactions_df[i,2]]
    }
  }
  
  names(reactions_df) <- c("source","target")
  
  reactions_df <- unique(reactions_df)
  
  ###
  
  row.names(reactions_df) <- c(1:length(reactions_df[,1]))
  
  if("Urea cycle" %in% pathway_to_keep)
  {
    reactions_df <- as.data.frame(rbind(reactions_df,network_supplements[["Urea cycle"]]))
  }
  if("Purine synthesis" %in% pathway_to_keep)
  {
    reactions_df <- as.data.frame(rbind(reactions_df,network_supplements[["Purine synthesis"]]))
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
  
  # compound_list <- c()
  # for( i in seq(1,length(compound_vec),10))
  # {
  #   compound_list <- c(compound_list, KEGGREST::keggGet(compound_vec[c(i:(i+9))]))
  # }
  
  bad_kegg_compounds <- list()
  for (compound in compound_list)
  {
    
    brite <- compound$BRITE
    if(!is.null(compound$ATOM))
    {
      brite <- gsub("[ ]+","",brite)
      if("Cofactors" %in% brite | "Nucleotides" %in% brite | "CO2;" %in% compound$NAME | as.numeric(compound$ATOM[1]) <= 3 | "ITP;" %in% compound$NAME | "IDP;" %in% compound$NAME | "NADH;" %in% compound$NAME | "NADPH;" %in% compound$NAME)
      {
        bad_kegg_compounds[[compound$ENTRY[1]]] <- compound
      }
    }
    # } else
    # {
    #   bad_kegg_compounds[[compound$ENTRY[1]]] <- compound
    # }
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
compress_transporters <- function(sub_network_nocofact)
{
  test_1 <- paste(gsub("_[cxrnme]$","",sub_network_nocofact$reaction_network$source), 
                  gsub("_[cxrnme]$","",sub_network_nocofact$reaction_network$target),
                  sep = "_")
  
  test_2 <- paste(gsub("_[cxrnme]$","",sub_network_nocofact$reaction_network$target), 
                  gsub("_[cxrnme]$","",sub_network_nocofact$reaction_network$source),
                  sep = "_")
  
  test_1[test_1 %in% test_2]
  
  transporters <- unique(c(sub_network_nocofact$reaction_network[test_1 %in% test_2,1], sub_network_nocofact$reaction_network[test_1 %in% test_2,2]))
  transporters <- transporters[!grepl("_[cxrnme]$",transporters)]
  
  for(i in 1:2)
  {
    sub_network_nocofact$reaction_network[,i] <- sapply(sub_network_nocofact$reaction_network[,i], function(x, transporters){
      if(x %in% transporters)
      {
        return("transporter")
      } else
      {
        return(x)
      }
    }, transporters = transporters)
  }
  
  sub_network_nocofact$reaction_network$edgeId <- paste(sub_network_nocofact$reaction_network$source, sub_network_nocofact$reaction_network$target, sep = "_")
  sub_network_nocofact$reaction_network <- sub_network_nocofact$reaction_network[!duplicated(sub_network_nocofact$reaction_network$edgeId),]
  
  edge_transporter <- sub_network_nocofact$reaction_network$edgeId
  edge_transporter <- edge_transporter[grepl("transporter",edge_transporter)]
  
  groups <- rep(0,length(edge_transporter))
  for(i in 1:length(edge_transporter))
  {
    if(grepl("^cpd:",edge_transporter[i]) | grepl("_[cxrnme]_",edge_transporter[i]))
    {
      metab <- gsub("_.*","",edge_transporter[i])
      for(j in 1:length(edge_transporter))
      {
        if(grepl("^transporter",edge_transporter[j]))
        {
          if(grepl(metab, edge_transporter[j]) & gsub("_transporter","",edge_transporter[i]) != gsub("transporter_","",edge_transporter[j]))
          {
            groups[i] <- i
            groups[j] <- i
          }
        }
      }
    }
  }
  
  names(groups) <- edge_transporter
  
  network <- sub_network_nocofact$reaction_network
  for(i in 1:2)
  {
    for(j in 1:length(network[,i]))
    {
      if(network[j,i] == "transporter")
      {
        print(network[j,3])
        network[j,i] <- paste(network[j,i], groups[network[j,3]], sep = ">")
      }
    }  
  }
  
  network <- network[!grepl("^SLC",network$edgeId) & !grepl("_SLC",network$edgeId),]
  
  sub_network_nocofact$reaction_network <- network[,-3]
  
  sub_network_nocofact$attributes <- sub_network_nocofact$attributes[sub_network_nocofact$attributes[,1] %in% network[,1] | sub_network_nocofact$attributes[,1] %in% network[,2],]
  
  new_transporters <- unique(c(network[,1],network[,2]))
  new_transporters <- new_transporters[grepl("transporter",new_transporters)]
  
  new_transporters <- as.data.frame(cbind(new_transporters,rep("transporter",length(new_transporters))))
  names(new_transporters) <- c("V1","V2")
  
  sub_network_nocofact$attributes <- as.data.frame(rbind(sub_network_nocofact$attributes, new_transporters))
  
  return(sub_network_nocofact)
}