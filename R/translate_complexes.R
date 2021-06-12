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
library(sjmisc)   


#' \code{translate_complexes}
#'
#' Function to separate the complexes into a vector of enzyme strings, map the
#' identifiers to their enzyme name using recon2_redhuman, and give back the full 
#' source/target network with enzymes, metabolites and translated complexes.
#'
#' @param Network source/target network containing enzymes, complexes, metabolites
#' @return source/target network containing enzymes, translated complexes, metabolites

translate_complexes <- function(Network){
  
  ##Add new column to network data frame for later comparisons
  Network$row_names <- 1:nrow(Network) 
  Network <- Network[c("row_names", "source", "target")]
  
  ##Select all strings that start with 3 digits "^\d{3}" (=complexes)
  network_source <- Network[grep("^\\d{3}", Network$source), ]
  network_target <- Network[grep("^\\d{3}", Network$target), ]  
  
  ##Delete all rows with complexes from network dataframe
  network_minus <- anti_join(Network, network_source, by= "row_names")
  network_minus <- anti_join(network_minus, network_target, by= "row_names")  
  
  ##Split complexes to separated enzymes using underscore as splitting condition
  split_complexes <- function(Column){
    for(i in 1:length(Column)){
      Column[i] <- strsplit(as.character(Column[i]), split = "_")
    }
    return(Column)
  }
  network_source$source <- split_complexes(network_source$source)
  network_target$target <- split_complexes(network_target$target)
  
  ##Combine source and target dataframes
  network_complex <- rbind(network_source, network_target)
  
  ##Use recon2_redhuman$gene_mapping and map enzyme names to identifiers
  mapping_vec <- recon2_redhuman$gene_mapping$name
  names(mapping_vec) <- recon2_redhuman$gene_mapping$X1
  
  for (i in 2:ncol(network_complex))      #iterate through columns 2 and 3 of df
  {
    for (j in 1:length(network_complex[,2]))      #j iterates through rows of df
    {
      str <- ""
      for (entry in network_complex[j,i][[1]])
      {
        if(entry %in% names(mapping_vec))
        {
          str <- paste(str, as.character(mapping_vec[[entry]]),sep="_")
        } else {
          str <- paste(str, as.character(entry),sep="_")
        }
        
        ##Take the ">" into account
        if(str_contains(entry,">"))
        {
          str <- gsub(entry, ">", str)
          str <- substr(str, 1, nchar(str)-3)
          pair <- strsplit(entry, split=">")
          
          tmp <- paste(mapping_vec[[pair[[1]][1]]],">",pair[[1]][2],sep="")
          str <- paste(str, tmp,sep="_")
        }
      }
      network_complex[j,i][[1]] <- substring(str,2)
    }
  }
  
  ##Add rows with translated complexes to network from beginning & order rows
  network_translated <- rbind(network_minus, network_complex)
  network_translated <- network_translated[order(network_translated[,1]), ]
  network_translated$row_names <- NULL
 
  ##To have the rows of the new order listed sequentially
  rownames(network_translated) <- 1:nrow(network_translated)
  
  ##Convert all columns to characters
  network_translated <- network_translated %>%
                              mutate(across(everything(), as.character))

  return(network_translated)
}