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


#'\code{translate_complexes}
#'
#'Function to separate the complexes into a vector of enzyme strings, map the
#'identifiers to their enzyme name using recon2_redhuman, and give back the full 
#'source/target network with enzymes, metabolites and translated complexes.
#'
#'@param  cellNetwork source/target network containing enzymes, complexes, metabolites
#'@return source/target network containing enzymes, translated complexes, metabolites

translate_complexes <- function(cellNetwork){
  
  
  # add new column to cell_network df for later comparisons
  cell_network$row_names <- 1:nrow(cell_network) 
  cell_network <- cell_network[c("row_names", "source", "target")]
  
  
  # select all strings that start with 3 digits "^\d{3}" (=complexes)
  # and delete those rows from cell_network dataframe
  cell_network_source <- cell_network[grep("^\\d{3}", cell_network$source), ]
  cell_network_minus <- anti_join(cell_network, cell_network_source, by= "row_names")
  # split enzymes in complexes to separate strings using underscore as splitting condition
  for(i in 1:length(cell_network_source$source)){
    cell_network_source$source[i] <- strsplit(as.character(cell_network_source$source[i]), split = "_")
  }
  
  cell_network_target <- cell_network[grep("^\\d{3}", cell_network$target), ]
  cell_network_minus <- anti_join(cell_network_minus, cell_network_target, by= "row_names")
  for(i in 1:length(cell_network_target$target)){
    cell_network_target$target[i] <- strsplit(as.character(cell_network_target$target[i]), split = "_")
  }
  
  # combine both dataframes
  cell_network_complex <- rbind(cell_network_source, cell_network_target)
  
  # 3. use recon2_redhuman$gene_mapping and map complexes to identifiers
  mapping_vec <- recon2_redhuman$gene_mapping$name
  names(mapping_vec) <- recon2_redhuman$gene_mapping$X1
  
  a <- cell_network_complex
  
  for (i in 2:ncol(a))   # columns 2 and 3 of df
  {
    for (j in 1:length(a[,2])) # j goes through rows of df
    {
      str <- ""
      for (entry in a[j,i][[1]])
      {
        if(entry %in% names(mapping_vec))
        {
          str <- paste(str, as.character(mapping_vec[[entry]]),sep="_")
        } else {
          str <- paste(str, as.character(entry),sep="_")
        }
        
        # Take the > into account
        if(str_contains(entry,">"))
        {
          print(entry)
          str <- gsub(entry, ">", str)
          str <- substr(str, 1, nchar(str)-3)
          pair <- strsplit(entry, split=">")
          
          value1 <- mapping_vec[[pair[[1]][1]]]
          value2 <- pair[[1]][2]
          
          tmp <- paste(value1,">",value2,sep="")
          str <- paste(str, tmp,sep="_")
        }
      }
      a[j,i][[1]] <- substring(str,2)
    }
  }
  
  # add rows with translated complexes to cell_network_minus
  cell_network_translated <- rbind(cell_network_minus, a)
  # order rows by column "row_names"
  cell_network_translated <- cell_network_translated[order(cell_network_translated[,1]), ]
  # delete column "row_names"
  cell_network_translated$row_names <- NULL
  # convert all columns to characters in dataframe
  cell_network_translated <- cell_network_translated %>%
                              mutate(across(everything(), as.character))
  
  return(cell_network_translated)
}