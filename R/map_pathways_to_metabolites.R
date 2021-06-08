############# Supporting function for translating enzyme complexes #############

#Copyright (C) 2021  Caroline Lohoff, Aurelien Dugourd
#Contact : aurelien.dugourd@bioquant.uni-heidelberg.de

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(ocean)
library(dplyr)
library(magrittr)
library(sjmisc)   # function "str_contains()"


#'\code{map_pathways_to_metabolites}
#'
#'Function to map pathways to enzymes and then to metabolites.
#'
#'@param  regulonsDf 
#'@return data frame containing metabolites and their pathways

map_pathways_to_metabolites <- function(regulonsDf){

  metab_complexes_df <- regulonsDf
  metab_complexes_df$weight <- NULL
  names(metab_complexes_df) <- c("enzymes", "metabolites")

  metab_df <- metab_complexes_df


  ##Split complexes: every enzyme gets a new row in df
  for (i in 1:nrow(metab_df))          # loop through all rows in "enzymes" column
  {
    if(str_contains(metab_df$enzymes[i],"_") & !str_contains(metab_df$enzymes[i],"_reverse"))
    {
      split_str <- as.vector(unlist(strsplit(as.character(metab_df$enzymes[i]), split = "_")))
    
      # create new row for every element (from the 2nd element on) in split_str
      for (j in 2:length(split_str))
      {
        new_row <- c(split_str[j], metab_df$metabolites[i])
        metab_df <- rbind(metab_df, new_row)
      }
      metab_df$enzymes[i] <- split_str[1]
    }
  }

  # reorder rows by two conditions: 1. column "enzymes", 2. column "metabolites"
  metab_df <- metab_df[order(metab_df$enzymes, metab_df$metabolites), ]
  # keep only unique rows (remove repetitions)
  metab_df_unique = metab_df %>% distinct()

  ### Is the information ">1234" and "_reverse" important when pathways are mapped to metabolic enzymes? 
  ### --> NO, because when we only have metabolites & pathways, we don't need the information!
  # MAKE LOOP FASTER

  metab_df[ , "pathway"] <- NA      # add new empty column "pathway" to dataframe

  for (i in 1:nrow(metab_df))          # loop through all rows in "enzymes" column
  {
    print(metab_df$enzymes[i])
    
    if(str_contains(metab_df$enzymes[i],">"))
    {
      pair <- strsplit(metab_df$enzymes[i], split=">")
    
      for (j in 1:nrow(kegg_pathways))
      {
        metab_df$pathway[i][pair[[1]][1] == kegg_pathways$gene[j]] <- (metab_df$pathway[i][pair[[1]][1] == kegg_pathways$gene[j]] <- kegg_pathways$term[j])
      }
      metab_df$enzymes[i] <- paste(pair[[1]][1],">",pair[[1]][2],sep="")
    }
    if(str_contains(metab_df$enzymes[i],"_reverse"))
    {
      metab_df$enzymes[i] <- gsub("_reverse","",metab_df$enzymes[i])
      metab_df$pathway[i][metab_df$enzymes[i] == kegg_pathways$gene[j]] <- (metab_df$pathway[i][metab_df$enzymes[i] == kegg_pathways$gene[j]] <- kegg_pathways$term[j])
      metab_df$enzymes[i] <- paste(metab_df$enzymes[i],"_reverse",sep="")
    }
  
    for (j in 1:nrow(kegg_pathways))
    {
      if (metab_df$enzymes[i] %in% kegg_pathways$gene[j])
      {
        metab_df$pathway[i][metab_df$enzymes[i] == kegg_pathways$gene[j]] <- (metab_df$pathway[i][metab_df$enzymes[i] == kegg_pathways$gene[j]] <- kegg_pathways$term[j])
      }
    }
    print(metab_df$pathway[i])
  }

  ## FASTER ALTERNATIVE (without considering ">1234" and "_reverse")

  ## delete ">" and digits behind
  #metab_df$enzymes <- gsub(">.*","",metab_df$enzymes)
  #metab_df$enzymes <- gsub("_reverse","",metab_df$enzymes) # delete "_reverse"

  # add new column pathways by merging metab df with kegg pathways (by enzymes)
  #metab_pathway_df = merge(metab_df[, c("enzymes", "metabolites")], 
  #                             kegg_pathways[, c("gene", "term")], # select all columns to be merged in df
  #                             by.x = "enzymes",         # select columns in each df by which merging should be done
  #                             by.y = "gene",
  #                             all.x = TRUE)      # keep all rows from 1st df even if there are no data in 2nd df
  # rename pathway column
  #colnames(metab_pathway_df)[colnames(metab_pathway_df) == "term"] <- "pathway"


  # change order of columns
  metab_pathway_df <- metab_pathway_df[c("metabolites",
                                       "metabolic_enzymes", "pathway")]
  # reorder rows by column "metabolites"
  metab_pathway_df <- metab_pathway_df[order(metab_pathway_df$metabolites), ]

  # keep only unique rows
  metab_pathway_dfu <- distinct(metab_pathway_df, .keep_all = TRUE)
  #pathways_nodes_unique <- pathways_nodes_df[!duplicated(pathways_nodes_df$Nodes),]


  # add new column deregulation_metabolite by merging metab df with t_table (by metabolites)
  metab_pathway_reg_df = merge(metab_pathway_df[, c("metabolites", "enzymes", "pathway")], 
                             translated_results$t_table[, c("KEGG", "tumorVsHealthy")], # select all columns to be merged in df
                             by.x = "metabolites",         # select columns in each df by which merging should be done
                             by.y = "KEGG",
                             all.x = TRUE) 
  # rename deregulation column
  colnames(metab_pathway_reg_df)[colnames(metab_pathway_reg_df) == "tumorVsHealthy"] <- "deregulation_metabolite"

  # keep only unique rows (remove repetitions)
  metab_pathway_reg_df = metab_pathway_reg_df%>%group_by_all

  # dataframe containing only metabolites and their pathways
  metabolites_pathways <- metab_pathway_reg_df[, c("metabolites", "pathway")]

  # save output files as RData
  saveRDS(metab_pathway_reg_df,"./results/metabolites_enzymes_pathways.RData")
  saveRDS(metabolites_pathways,"./results/metabolites_pathways.RData")

  return(metabolites_pathways)
}