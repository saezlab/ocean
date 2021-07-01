#' target_set_from_forest_2
#'
#' this function converts a 'forest' representation of a metabolic network into a metabolic-enzyme set
#'
#' @param forest  a forest representation of a metabolic network. can be used with the tree_without_cofactor
#' @param measured_species character vector, the names of the metabolites that are measured. 
#' Must be coherents with identifers in forest
#' @param penalty numeric, the penalty, between 0.1 and 0.9 (0 and 1 work as well but are special cases and shouldn't be used)
#' @param verbose boolean, to display progress
#' @return a reaction network in the form of a three column dataframe with enzme, metabolite targets and weights.
target_set_from_forest_2 <- function(forest, measured_species, penalty = 0.5, verbose = T)
{
  target_set <- list()
  i <- 1
  l <- 1
  for(tree in forest)
  {
    if(verbose)
    {
      print(paste0(i, "/", length(names(forest)), " ", "pen = ",penalty))
    }
    direction_sign <- -1
    for (direction in tree)
    {
      j <- 1
      for (layer in direction)
      {
        if(length(layer) > 0)
        {
          layer_set <- as.data.frame(matrix(NA,length(layer), 3))
          names(layer_set) <- c("set","targets","weight")

          layer_set$set <- names(forest)[i]
          layer_set$targets <- layer
          layer_set$weight <- j*direction_sign

          j <- j*penalty

          layer_set <- layer_set[layer_set$targets %in% measured_species,]

          if(length(layer_set[,1] > 0))
          {
            target_set[[l]] <- layer_set
            l <- l+1
          }
        }
      }
      direction_sign <- direction_sign*-1
    }
    i <- i+1
  }
  target_set <- as.data.frame(do.call(rbind,target_set))
  return(target_set)
}

#' prepare_metabolite_set
#'
#' This functions create a series of metabolic-enzyme set collections over a range of penalties
#' (see tutorial)
#'
#' @param penalty_range  numeric vector, between minimum 1 and maximum 9
#' @param forest list, a forest representation of a metabolic network. 
#' @param measured_metabolites character vector, names of measured metabolites
#' @return a series of metabolic-enzyme set collections
#' @export
prepare_metabolite_set <- function(penalty_range, forest, measured_metabolites)
{
  reaction_set_list <- list()
  for(i in penalty_range)
  {
    print(i)
    pen <- i/10
    reaction_set_list[[i]] <- target_set_from_forest_2(forest, measured_metabolites, penalty = pen)
  }
  return(reaction_set_list)
}

#' condense_metabolite_set
#'
#' this function condense the metabolic-enzyme set collections to summarise ignore specific types of
#' information such as compartments of relative positions (sign)
#'
#' @param reaction_set_list  list, metabolic-enzyme set collections as returned by the prepare_metabolite_set function
#' @param condense_sign boolean, summarise upstream and downstream position of metabolites into a signle relative position
#' @param condense_compartments boolean, wether the compartment information of hte metabolic network sohuld be ignored
#' @param ignore_sign boolean, wether the upstream/downstream position of metabolic should be ignored 
#' @return a metabolic-enzyme set collections
#' @export
#' @import dplyr
#' @import magrittr
condense_metabolite_set <- function(reaction_set_list, condense_sign = T, condense_compartments = F, ignore_sign = F)
{
  reaction_set_list_merged <- list()
  for(i in 1:length(reaction_set_list))
  {
    if(!is.null(reaction_set_list[[i]]))
    {
      reaction_set <- reaction_set_list[[i]]

      #to merge sign
      if(condense_sign)
      {
        reaction_set <- reaction_set %>% group_by(set,targets) %>% summarise_each(funs(sum(., na.rm = TRUE)))
        reaction_set <- as.data.frame(reaction_set)
      }

      #to merge compartment
      if(condense_compartments)
      {
        reaction_set[,2] <- gsub("_[a-z]$","",reaction_set[,2])
        reaction_set <- reaction_set %>% group_by(set,targets) %>% summarise_each(funs(mean(., na.rm = TRUE)))
        reaction_set <- as.data.frame(reaction_set)
      }


      #to ignore upstream/downstream
      if(ignore_sign)
      {
        reaction_set$weight <- abs(reaction_set$weight)
      }

      reaction_set_list_merged[[i]] <- reaction_set
    }
    else
    {
      reaction_set_list_merged[[i]] <- reaction_set_list[[i]]
    }
  }
  return(reaction_set_list_merged)
}

#' prepare_regulon_df
#'
#' extract and cleanup the regulon collection correpsonding to the desired penalty
#'
#' @param reaction_set_list_merged  a metabolic-enzyme set collections, 
#' as returned by the condense_metabolite_set function
#' @param penalty numeric, the index of the reaction_set_list_merged corresponding to
#' the desired penalty
#' @param filter_imbalance a vector of two numbers, representing the lower and upper bounds to filter 
#' imbalanced reaction trees. for example c(0.33,0.67) will filter out any reaction were the upstream branch is 
#' longer than 2 time the length of the downstream branch, and vice versa
#' @return a regulon dataframe with three columns correposning to the enzme, the target metabolite and the weight
#' @export
prepare_regulon_df <- function(reaction_set_list_merged, penalty, filter_imbalance = c(0,1))
{
  n <- length(unique(reaction_set_list_merged[[penalty]][,1]))

  regulons_df <- reaction_set_list_merged[[penalty]]

  regulons_df <- regulons_df[regulons_df[,1] %in% unique(regulons_df[,1])[1:n],]
  regulons_df <- unique(regulons_df)
  
  test <- sapply(unique(regulons_df$set), function(x, regulons_df){
    sub_reg <- regulons_df[regulons_df$set == x,]
    prop <- sum(sub_reg$weight < 0) / length(sub_reg$weight)
    if(prop >= filter_imbalance[1] & prop <= filter_imbalance[2])
    {
      return(x)
    } else
    {
      return(NA)
    }
  }, regulons_df = regulons_df)
  
  regulons_df <- regulons_df[regulons_df$set %in% test,]

  return(regulons_df)
}

#' metactivity
#'
#' Use a mean enrichement statistic with shuffling-based normalisation (preserving metabolic compartment information)
#' to compute the metabolic enzme normalised enrichment score, based on its upstream and downstream metabolites
#'
#' @param metabolomic_t_table  ipsum...
#' @param regulons_df ipsum...
#' @param compartment_pattern ipsum...
#' @param k ipsum...
#' @return ipsum...
#' @export
#' @importFrom reshape2 dcast
metactivity <- function(metabolomic_t_table,
                        regulons_df,
                        compartment_pattern = "",
                        k = 1000,
                        seed = 1337)
{
  set.seed(seed)
  enzymes_nes_list <- list()
  enzymes_es_list <- list()
  for(i in 2:length(metabolomic_t_table[1,]))
  {

    if(compartment_pattern != "")
    {
      t_table <- metabolomic_t_table
      t_table$unique_metab <- gsub(compartment_pattern, "", metabolomic_t_table[,1])
      t_table_reduced <- t_table[,c(length(t_table[1,]),2:(length(t_table[1,])-1))]
      t_table_reduced <- unique(t_table_reduced)

      t_null <- as.data.frame(replicate(k, sample(t_table_reduced[,i], length(t_table_reduced[,i]))))
      t_null <- as.data.frame(cbind(t_table_reduced[,1],t_null))
      names(t_null)[1] <- names(t_table_reduced)[1]

      t_table_with_null <- merge(t_table[,c(1,i,length(t_table[1,]))],
                                 t_null,
                                 by = names(t_table_reduced)[1])
      t_table_with_null <- t_table_with_null[,-1]
      row.names(t_table_with_null) <- t_table_with_null[,1]
      t_table_with_null <- t_table_with_null[,-1]
    } else
    {
      t_table <- metabolomic_t_table
      t_table$unique_metab <- metabolomic_t_table[,1]
      t_table_reduced <- t_table[,c(length(t_table[1,]),2:(length(t_table[1,])-1))]
      t_table_reduced <- unique(t_table_reduced)
      
      t_null <- as.data.frame(replicate(k, sample(t_table_reduced[,i], length(t_table_reduced[,i]))))
      t_null <- as.data.frame(cbind(t_table_reduced[,1],t_null))
      names(t_null)[1] <- names(t_table_reduced)[1]
      
      t_table_with_null <- merge(t_table[,c(1,i,length(t_table[1,]))],
                                 t_null,
                                 by = names(t_table_reduced)[1])
      t_table_with_null <- t_table_with_null[,-1]
      row.names(t_table_with_null) <- t_table_with_null[,1]
      t_table_with_null <- t_table_with_null[,-1]
    }

    regulons_df <- regulons_df[regulons_df[,2] %in% row.names(t_table_with_null),]

    regulons_mat <- reshape2::dcast(regulons_df, targets~set, value.var = "weight")
    row.names(regulons_mat) <- regulons_mat[,1]
    regulons_mat <- regulons_mat[,-1]

    t_table_with_null <- t_table_with_null[row.names(t_table_with_null) %in% row.names(regulons_mat),]

    regulons_mat <- regulons_mat[row.names(t_table_with_null),]
    regulons_mat[is.na(regulons_mat)] <- 0


    metabolites <- row.names(t_table_with_null)
    enzymes <- names(regulons_mat)
    t_table_with_null <- t(t_table_with_null)

    enzyme_ES <- as.matrix(t_table_with_null) %*% as.matrix(regulons_mat)

    enzymes_es_list[[i]] <- enzyme_ES[1,]

    enzyme_NES <- enzyme_ES[1,]
    null_means <- colMeans(enzyme_ES[-1,])
    null_sds <- apply(enzyme_ES[-1,],2,sd)

    enzyme_NES <- (enzyme_NES - null_means) / null_sds

    enzymes_nes_list[[i]] <- enzyme_NES
  }

  enzymes_nes_df <- as.data.frame(do.call(cbind,enzymes_nes_list))
  enzymes_nes_df$enzyme <- row.names(enzymes_nes_df)
  enzymes_nes_df <- enzymes_nes_df[,c(length(enzymes_nes_df[1,]),1:(length(enzymes_nes_df[1,])-1))]
  names(enzymes_nes_df) <- names(metabolomic_t_table)

  enzymes_es_df <- as.data.frame(do.call(cbind,enzymes_es_list))
  enzymes_es_df$enzyme <- row.names( enzymes_es_df)
  enzymes_es_df <-  enzymes_es_df[,c(length( enzymes_es_df[1,]),1:(length( enzymes_es_df[1,])-1))]
  names( enzymes_es_df) <- names(metabolomic_t_table)

  return(list("ES" = enzymes_es_df, "NES" = enzymes_nes_df))
}
