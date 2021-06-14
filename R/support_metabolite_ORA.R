############## Supporting functions for ORA based on metabolites ###############

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
library(sjmisc)   
library(ggplot2)


#' \code{rearrange_dataframe}
#'
#' Function to separate the complexes into a vector of enzyme strings, map the
#' identifiers to their enzyme name using recon2_redhuman, and give back the full 
#' source/target network with enzymes, metabolites and translated complexes.
#'
#' @param Network source/target network containing enzymes and metabolites
#' @return data frame with metabolites and their target/source enzymes

rearrange_dataframe <- function(network){
  
  new_cols = c("enzymes", "metabolites")
  network[new_cols] <- NA
  metabolites_df <- network[c("enzymes", "metabolites", "source", "target")]
  
  ##Copy entries from source and target column to metabolite or enzyme column
  
  for (i in 3:ncol(metabolites_df)) #iterate through source and target column
  {
    for (j in 1:nrow(metabolites_df)) #iterate through rows
    {
      if (str_contains(metabolites_df[j,i],"cpd:")){
        metabolites_df$metabolites[j] <- metabolites_df[j,i]
      } else
      {
        metabolites_df$enzymes[j] <- metabolites_df[j,i]
      }
    }
  }
  
  metabolites_df <- metabolites_df[, -c(3:4)]  #delete columns source and target
  metabolites_df <- na.omit(metabolites_df)    #remove rows containing NA
  
  return(metabolites_df)
}



#'\code{map_pathways_to_metabolites}
#'
#'Function to map pathways to enzymes and then to metabolites.
#'
#'@param  metab_df data frame with metabolites and their target/source enzymes 
#'@return data frame containing metabolites and their pathways

map_pathways_to_metabolites <- function(metab_df){
  
  metab_df$enzymes <- gsub(">.*","",metab_df$enzymes) #delete ">" and digits behind
  metab_df$enzymes <- gsub("_reverse","",metab_df$enzymes) #delete "_reverse"
  metab_df <- distinct(metab_df)  #keep only unique rows
  
  ##Split complexes: every enzyme gets a new row in df
  for (i in 1:nrow(metab_df))    #iterate through all rows in "enzymes" column
  {
    if(str_contains(metab_df$enzymes[i],"_"))
    {
      split_str <- as.vector(unlist(strsplit(as.character(metab_df$enzymes[i]),
                                             split = "_")))
      
      #create new row for every element in split_str (from the 2nd element on) 
      for (j in 2:length(split_str))
      {
        new_row <- c(split_str[j], metab_df$metabolites[i])
        metab_df <- rbind(metab_df, new_row)
      }
      metab_df$enzymes[i] <- split_str[1]
    }
  }
  
  ##Reorder rows by two conditions: 1. column "enzymes", 2. column "metabolites"
  metab_df <- metab_df[order(metab_df$enzymes, metab_df$metabolites), ]
  ##Keep only unique rows (if a complex consisted of two or more similar enzymes)
  metab_df <- distinct(metab_df)
  

  ##Add new column "pathways" by merging data frame with kegg pathways (by enzymes)
  metab_pathway_df = merge(metab_df[, c("enzymes", "metabolites")], 
                           kegg_pathways[, c("gene", "term")], 
                           by.x = "enzymes",         
                           by.y = "gene",
                           all.x = TRUE)      

  colnames(metab_pathway_df)[colnames(metab_pathway_df) == "term"] <- "pathway"
  
  metab_pathway_df$enzymes <- NULL  # remove column "enzymes"
  
  ##Reorder rows by column "metabolites"
  metab_pathway_df <- metab_pathway_df[order(metab_pathway_df$metabolites), ]
  metab_pathway_df <- na.omit(metab_pathway_df)   #remove rows containing NA
  metab_pathway_df <- distinct(metab_pathway_df)  #keep only unique rows
  
  return(metab_pathway_df)
}



#'\code{barplot_pathways}
#'
#'Function to plot the pathways which are significantly overexpressed
#'
#'@param sigPathwaysDf  Data frame containing all significant pathways as a result of the ORA
#'@return Bar plot of the most significant pathways by -log10(p-value)

barplot_pathways <- function(sigPathwaysDf){
  
  ##Filtering out top pathways (max. 20) by Adjusted p-value
  sigPathwaysDf <- sigPathwaysDf[order(sigPathwaysDf$'Adjusted p-value'),]
  top_hallmark <- sigPathwaysDf[1:20, c(7,1,2)]
  top_hallmark <- na.omit(top_hallmark)   #remove rows containing NA
  top_hallmark <- top_hallmark[order(top_hallmark$'p-value', decreasing = TRUE),]
  
  ##Write pathways with capital letter and without underscore
  top_hallmark$pathway <- tolower(top_hallmark$pathway)
  top_hallmark$pathway <- gsub("(^|\\p{P})(.)",
                               "\\1\\U\\2",    
                               top_hallmark$pathway, perl = TRUE)
  top_hallmark$pathway <- gsub("_"," ",top_hallmark$pathway)
  
  top_hallmark$`p-value` <- -log10(top_hallmark$`p-value`)
  top_hallmark$pathway <- factor(top_hallmark$pathway, levels = top_hallmark$pathway)
  names(top_hallmark)[3] <- "-log10(p-value)"
  
  ##Create bar plot
  plot <- ggplot(top_hallmark, aes(x = pathway, y = `-log10(p-value)`,
                                   fill = `-log10(p-value)`)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    theme_minimal() + 
    ggtitle("ORA Metabolites Top Pathways") +
    scale_fill_gradient(low="grey", high="darkred")
  
#  ggsave("top_pathways.png", plot = plot,
#         path = "results/", scale = 1, dpi = 300, limitsize = TRUE)
  
  return(plot)
}
