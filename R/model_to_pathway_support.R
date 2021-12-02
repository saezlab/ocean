build_grps <- function(rxns, rules, minmax, pathways)
{
  gprs <- as.data.frame(cbind(rxns,rules))
  
  names(gprs) <- c("Reactions","Genes")
  
  gprs$Genes <- gsub("[()]","", gprs$Genes)
  
  gprs$Genes <- gsub(" or ",";",gprs$Genes)
  gprs$Genes <- gsub(" and ","_",gprs$Genes)
  gprs$Genes <- gsub("[.][0-9]","",gprs$Genes)
  
  gprs$direction <- minmax$direction
  
  gprs$pathway <- pathways
  
  return(gprs)
}

build_reactions_list <- function(gprs, stochio, metab, minmax)
{
  reactions_list <- list()
  for(i in 1:length(gprs[,1]))
  {
    reactions_list[[i]] <- list()
    reactants <- which(stochio[,i] <= -1)
    
    if(!identical(reactants, integer(0)))
    {
      reactants <- metab[reactants,]
      reactions_list[[i]]$reactants <- reactants
    }
    else
    {
      reactions_list[[i]]$reactants <- NA
    }
    
    products <- which(stochio[,i] >= 1)
    
    if(!identical(products, integer(0)))
    {
      products <- metab[products,]
      reactions_list[[i]]$products<- products
    }
    else
    {
      reactions_list[[i]]$products <- NA
    }
    reactions_list[[i]]$direction <- minmax[i,3]
  }
  names(reactions_list) <- gprs$Reactions
  
  for(i in 1:length(reactions_list))
  {
    reaction <- names(reactions_list)[i]
    
    reactions_list[[i]]$genes <- gprs[gprs$Reactions == reaction,2]
  }
  
  return(reactions_list)
}

build_reactions_df <- function(reactions_list, mapping)
{
  mapping_vec <- mapping$name
  names(mapping_vec) <- mapping$X1
  
  reactions_df <- as.data.frame(matrix(NA,0,2)) 
  
  genes_memory <- c()
  z <- 1
  m <- 1
  for(reaction in reactions_list)
  {
    #First, get the genes, metabolite reactant and products that are involved in the current reaction
    if(!is.null(reaction$reactants) &!is.null(reaction$products))
    {
      if(!is.null(reaction) & !is.na(reaction$genes))
      {
        genes <- strsplit(reaction$genes,";")
        genes <- genes[[1]]
      }
      else
      {
        genes <- names(reactions_list)[m]
      }
      
      if(reaction$direction == 1)
      {
        reaction$reversible <- FALSE
        reactants <- reaction$reactants
        products <- reaction$products
      }
      else
      {
        if(reaction$direction != 0) #The reaction is operating in reversed mode
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
      
      for(gene in genes)
      {
        
        if(!(gene %in% genes_memory))
        {
          genes_memory <- c(genes_memory, gene)
          if(gene %in% names(mapping_vec))
          {
            gene <- mapping_vec[gene]
          }
        }
        else
        {
          if(gene %in% names(mapping_vec))
          {
            gene <- mapping_vec[gene]
          }
          gene <- paste(gene,z,sep = ">")
          # genes_memory <- c(genes_memory, gene)
          z <- z+1
        }
        
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
        
        reaction_df <- as.data.frame(rbind(reaction_df,new_reaction_df))
      }
      
      reactions_df <- as.data.frame(rbind(reactions_df, reaction_df))
    }
    m <- m+1
  }
  return(reactions_df)
}

# translate_network_gene_ids <- function(reactions_df, mapping)
# {
#   mapping_vec <- mapping$name
#   names(mapping_vec) <- mapping$X1
#   
#   for (i in 1:2)
#   {
#     for (j in 1:length(reactions_df[,1]))
#     {
#       if(grepl(">",reactions_df[j,i]))
#       {
#         affixe <- gsub(".+[>]","",reactions_df[j,i])
#         reactions_df[j,i] <- gsub("[>].*","",reactions_df[j,i])
#         has_affixe <- TRUE
#       }
#       else
#       {
#         if(grepl("<",reactions_df[j,i]))
#         {
#           affixe_lower <- gsub(".+[<]","",reactions_df[j,i])
#           reactions_df[j,i] <- gsub("[<].*","",reactions_df[j,i])
#           has_affixe_lower <- TRUE
#         }
#         else
#         {
#           has_affixe_lower <- FALSE
#         }
#         has_affixe <- FALSE
#       }
#       if(reactions_df[j,i] %in% names(mapping_vec))
#       {
#         if(has_affixe)
#         {
#           reactions_df[j,i] <- paste(mapping_vec[reactions_df[j,i]],affixe,sep = ">")
#         }
#         else
#         {
#           if(has_affixe_lower)
#           {
#             reactions_df[j,i] <- paste(mapping_vec[reactions_df[j,i]],affixe_lower,sep = "<")
#           }
#           else
#           {
#             reactions_df[j,i] <- mapping_vec[reactions_df[j,i]]
#           }
#         }
#       } else
#       {
#         if(has_affixe)
#         {
#           reactions_df[j,i] <- paste(reactions_df[j,i],affixe,sep = ">")
#         }
#         else
#         {
#           if(has_affixe_lower)
#           {
#             reactions_df[j,i] <- paste(reactions_df[j,i],affixe_lower,sep = "<")
#           }
#           else
#           {
#             reactions_df[j,i] <- reactions_df[j,i]
#           }
#         }
#       }
#       if(gsub("_reverse","",reactions_df[j,i]) %in% names(mapping_vec))
#       {
#         reactions_df[j,i] <- paste(mapping_vec[gsub("_reverse","",reactions_df[j,i])],"_reverse",sep = "")
#       }
#     }
#   }
#   return(reactions_df)
# }

translate_network_metab_ids <- function(reactions_df, mapping)
{
  mapping_named_vec <- mapping[,2]
  names(mapping_named_vec) <- mapping[,1]
  
  names(mapping_named_vec) <- gsub("__","_",names(mapping_named_vec))
  
  for (i in 1:length(reactions_df[,1]))
  {
    if (reactions_df[i,1] %in% names(mapping_named_vec))
    {
      reactions_df[i,1] <- mapping_named_vec[reactions_df[i,1]]
    }
    if (reactions_df[i,2] %in% names(mapping_named_vec))
    {
      reactions_df[i,2] <- mapping_named_vec[reactions_df[i,2]]
    }
  }
  
  return(reactions_df)
}