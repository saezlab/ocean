#' translate_results
#'
#' ipsum...
#'
#' @param regulons_df  ipsum...
#' @param t_table ipsum...
#' @param mapping_table ipsum...
#' @return ipsum...
#' @export
#' @importFrom stringr str_extract
#' @import dplyr
translate_results <- function(regulons_df, t_table, mapping_table, compress_compartments = F)
{
  mapping_vec <- mapping_table$metab
  names(mapping_vec) <- mapping_table$KEGG
  
  regulons_df_with_names <- regulons_df
  regulons_df_with_names$targets <- sapply(regulons_df_with_names$targets, function(x,mapping_vec){
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
    
  }, mapping_vec = mapping_vec, simplify = F, USE.NAMES = F)
  
  t_table_with_names <- t_table
  
  t_table_with_names$KEGG <- sapply(t_table_with_names$KEGG, function(x,mapping_vec){
    suffixe <- str_extract(x, "_.$")
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
    
  }, mapping_vec = mapping_vec, simplify = F, USE.NAMES = F)
  
  t_table_with_names$KEGG <- as.character(t_table_with_names$KEGG)
  regulons_df_with_names$targets <- as.character(regulons_df_with_names$targets)
  
  if(compress_compartments)
  {
    regulons_df_with_names$ID <- paste(regulons_df_with_names$set, gsub("_[a-z]$","",regulons_df_with_names$targets), sep = "___")
    regulons_df_with_names <- regulons_df_with_names[,-c(1,2)]
    
    regulons_df_with_names <- regulons_df_with_names %>% dplyr::group_by(ID) %>% dplyr::summarise_each(funs(mean(., na.rm = TRUE)))
    regulons_df_with_names <- as.data.frame(regulons_df_with_names)
    
    regulons_df_with_names$set <- gsub("___.*","",regulons_df_with_names$ID)
    regulons_df_with_names$targets <- gsub(".*___","",regulons_df_with_names$ID)
    regulons_df_with_names <- regulons_df_with_names[,c(3,4,2)]
  }
  
  return(list('t_table' = t_table_with_names, "regulons_df" = regulons_df_with_names))
}

#' plotMetaboliteContribution
#'
#' ipsum...
#'
#' @param enzyme  ipsum...
#' @param stat_df ipsum...
#' @param metabolite_sets ipsum...
#' @param contrast_index ipsum...
#' @param stat_name ipsum...
#' @param scaling_factor ipsum...
#' @param nLabels ipsum...
#' @return ipsum...
#' @export
#' @import ggplot2
#' @import ggrepel
plotMetaboliteContribution <- function(enzyme, stat_df, metabolite_sets, contrast_index, stat_name = "stat", scaling_factor = 1, nLabels = 10)
{
  
  metabolite_sets <- metabolite_sets[metabolite_sets$set == enzyme,]
  
  stat_df <- stat_df[stat_df[,1] %in% metabolite_sets$targets,]
  stat_df <- stat_df[,c(1,contrast_index+1)]
  
  
  stat_df <- unique(stat_df)
  
  names(stat_df) <- c("ID","C1")
  names(metabolite_sets) <- c("set","ID","weight")
  
  df <- merge(metabolite_sets,stat_df, by = "ID")
  
  df$contribution <- df$weight*df$C1*scaling_factor
  # df$contribution <- rescale(df$contribution, to = c(-3,3))
  df <- df[order(abs(df$contribution),decreasing = T),]
  
  save_labels <- df$ID  
  
  if(length(df$ID) > nLabels)
  {
    df$ID <- c(df[1:nLabels,'ID'],rep(NA,length(df$ID) - nLabels)) 
  }
  
  # df$ID <- ifelse(abs(df$contribution) >=  1,df$ID,"")
  # df$ID <- ifelse(abs(df$C1) >=  2,df$ID,"")
  
  ggscat <- ggplot(df, aes(x = weight, y = C1, label = ID)) + 
    geom_point(aes(fill = contribution),colour = "black",pch=21, size = abs(df$contribution)) +
    xlim(c(-1,1)) +
    geom_text_repel(point.padding = unit(max(abs(df$contribution)), "points")) +
    # scale_color_scico(palette = "vik", limits = c(-1, 1) * max(abs(df$contribution))) +
    scale_fill_gradient2(low="blue", high="red", midpoint = 0, mid = "white") +
    theme_minimal() + 
    geom_abline(intercept = 0, slope = 0, color = "grey") +
    geom_vline(xintercept = 0, color = "grey") +
    ylim(c(-max(abs(df$C1)),max(abs(df$C1)))) + 
    ggtitle(paste0(enzyme," metabolic consumption/production profile")) +
    labs(x = "consuption <==> production", y = stat_name)
  
  df <- df[order(df$contribution,decreasing = F),]
  df$running_sum_contribution <- cumsum(df$contribution)
  
  df$ID <- save_labels
  df$ID <- paste0(c(1:length(df$ID)), " ", df$ID)
  df$ID <- factor(df$ID, level = df$ID)
  df$group <- "A"
  df$label <- gsub("[0-9]+ ","",df$ID)
  
  ylim_bot <- ifelse(min(df$running_sum_contribution) < 0, min(df$running_sum_contribution)*1.25, 0 ) 
  ylim_bot <- ifelse(ylim_bot < min(df$contribution), ylim_bot, min(df$contribution) * 1.25)
  
  ylim_top <- ifelse(max(df$running_sum_contribution) > 0, max(df$running_sum_contribution)*1.25, 0 ) 
  ylim_top <- ifelse(ylim_top > max(df$contribution), ylim_top, max(df$contribution) * 1.25)
  
  ggcumSum <- ggplot(df, aes(x = ID, y = running_sum_contribution, group = group, label = label)) + 
    geom_point(aes(x = ID, y = contribution)) +
    geom_hline(yintercept = 0) +
    geom_line() + 
    theme_minimal() + 
    ylim(c(ylim_bot,ylim_top)) +
    scale_x_discrete(guide = guide_axis(n.dodge=8)) +
    theme(plot.margin = unit(c(2,2,2,2), "cm"))
  
  return(list("scatter" = ggscat, "cumsumPlot" = ggcumSum))
}

#' pathway_HM
#'
#' ipsum...
#'
#' @param mean_NES_df  ipsum...
#' @param pathway_name ipsum...
#' @param pathways ipsum...
#' @param sort_by ipsum...
#' @param manual_pathway ipsum...
#' @return ipsum...
#' @export
#' @importFrom reshape2 melt
pathway_HM <- function(mean_NES_df, pathway_name, pathways, sort_by = 1, manual_pathway = F)
{
  current_pathway <- pathways[pathways$term == pathway_name,"gene"]
  
  mean_NES_df_to_subset <- mean_NES_df
  
  if(manual_pathway == F)
  {
    mean_NES_df_to_subset$ID_short <- gsub(">.*","",mean_NES_df$KEGG)
    mean_NES_df_to_subset$ID_short <- gsub("_.*","",mean_NES_df_to_subset$ID_short)
  } else
  {
    mean_NES_df_to_subset$ID_short <- mean_NES_df_to_subset$KEGG
  }
  
  current_pathway_activities <- mean_NES_df[mean_NES_df_to_subset$ID_short %in% current_pathway,]
  
  row.names(current_pathway_activities) <- current_pathway_activities$KEGG
  current_pathway_activities <- current_pathway_activities[order(current_pathway_activities[,sort_by+1], decreasing = F),]
  
  
  current_pathway_act_melt <- melt(current_pathway_activities)
  current_pathway_act_melt$KEGG <- factor(current_pathway_act_melt$KEGG,levels = unique(current_pathway_act_melt$KEGG))
  names(current_pathway_act_melt) <- c("Enzyme","Contrast","NES")
  
  gp <- ggplot(current_pathway_act_melt, aes(x=Contrast, y = Enzyme, fill = NES)) + geom_tile() +
    scale_fill_gradient2(low="blue", high="red", midpoint = 0, mid = "white") + 
    theme_minimal() +
    labs(y = "Enzyme", x = "Contrast") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(gp)
}