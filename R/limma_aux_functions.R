#' checkInputs
#'
#' ipsum...
#'
#' @param measurments  ipsum...
#' @param targets ipsum...
#' @return ipsum...
checkInputs <- function(measurments, targets)
{
  if(class(measurments) != "data.frame")
  {
    error_message <- paste("The measurments argument should be a data.frame. It's currently a", paste(class(measurments), ".",sep = ""))
    return(list(FALSE, error_message))
  }
  else
  {
    if(dim(measurments)[1] == 0)
    {
      error_message <- "The measurments dataframe doesn't seem to contain any measurments..."
      return(list(FALSE, error_message))
    }
    else
    {
      if(dim(measurments)[2] == 0)
      {
        error_message <- "The measurments dataframe doesn't seem to contain any samples..."
        return(list(FALSE, error_message))
      }
      else
      {
        if(class(as.matrix(measurments)[,1]) != "numeric")
        {
          return(list(FALSE, "The measurments dataframe should contain only numerical values (or NAs)."))
        }
        else
        {
          if(class(targets) != "data.frame")
          {
            error_message <- paste("The targets argument should be a data.frame. It's currently a", paste(class(targets), ".",sep = ""))
            return(list(FALSE, error_message))
          }
          else
          {
            if(dim(targets)[2] < 2)
            {
              return(list(FALSE,"The targets dataframe should have at least two columns, sample names and conditions."))
            }
            else
            {
              if(dim(targets)[1] != dim(measurments)[2])
              {
                error_message <- paste("The targets dataframe should have as many samples (targets rows) as the measurements (measurments columns). Currently, the targets dataframe has", paste(dim(targets)[1], "samples and the measurements have", paste(dim(measurments)[2],"samples.")))
                return(list(FALSE, error_message))
              }
              else
              {

              }
            }
          }
        }
      }
    }
  }
  return(list(TRUE, "All seems to be in order..."))
}

#'\code{makeContrastsAlt}
#'
#'This function create a contrast matrix to be used by limma.
#'
#'@param targets A n*2 dataframe, where n is the number of samples. First column correspond to samples, second column correspond to conditions.
#'@param comparisons a list of numeric vectors. Each vector represent which condition should be conpared. Example :
#'c(2,-1) means that the first condition should be substracted from second condition. Vectors can be more than two element for complex contrasts.
#'@return a contrast matrix
makeContrastsAlt <- function(targets, comparisons)
{
  cont.matrix <- matrix(0,nrow = length(unique(targets$condition)), ncol = length(comparisons))
  i <- 1
  for (comparison in comparisons)
  {
    for (j in 1:length(comparison))
    {
      cont.matrix[abs(comparison[j]),i] <- cont.matrix[abs(comparison[j]),i]+(comparison[j]/abs(comparison[j]))
    }
    i <- i + 1
  }
  return(cont.matrix)
}

#'\code{runLimma}
#'
#'This function is a wrapper of limma made to facilitate the use of limma differential analysis
#'
#'@param measurments the measurment n*m dataframe (n is number of omic features, m is number of samples) where columns are ordered by conditions.
#'@param targets A n*2 dataframe, where n is the number of samples. First column correspond to samples, second column correspond to conditions.
#'@param comparisons a list of numeric vectors. Each vector represent which condition should be conpared. Example :
#'c(2,-1) means that the first condition should be substracted from second condition. Vectors can be more than two element for complex contrasts.
#'@param regress_out in case the user which to exclude possible confounding factors from the analysis, the user can provide additional columns in the targets dataframe.
#'then, the confounding factor can be regressed out by indicating the number of the column of the target dataframe describing it. Only one factor can be regressed out at the present time.
#'@return a list. First element is the limma model fitted with the contrast matrix, this is the usual output of limma. Second element is the contrast matrix that was used. third element is the fitted limma object without contrasts.
#' @export
#' @import limma
runLimma <- function(measurements, targets, comparisons = NULL, regress_out = NULL)
{
  input_check <- checkInputs(measurements, targets)
  if (input_check[[1]]) #input has correct format
  {
    if (!is.null(comparisons))
    {
      if (!is.null(regress_out))
      {
        for (regressor in regress_out)
        {
          measurements <- removeBatchEffect(measurements, targets[,regressor])
        }
      }

      cont.matrix <- makeContrastsAlt(targets, comparisons)

      cont.matrix <- as.data.frame(cont.matrix)
      row.names(cont.matrix) <- unique(targets$condition)
      cont.matrix <- as.matrix(cont.matrix)

      fcond <- factor(targets$condition, levels = unique(targets$condition))

      design <- model.matrix(~0+fcond)
      design <- as.data.frame(design)
      names(design) <- unique(targets$condition)
      design <- as.matrix(design)

      print(cont.matrix)

      fit <- limma::lmFit(measurements, design)
      fit2 <- limma::contrasts.fit(fit, cont.matrix)
      fit2 <- limma::eBayes(fit2)

      return(list(fit2, cont.matrix, fit))
    }
  }
  else
  {
    print(input_check[[2]])
    return(input_check[[1]])
  }
}

#' ttop_list_to_t_table
#'
#' This function converts a list of limma top table results into a single 
#' dataframe of t-values. The list of top table results can be obtained using 
#' the limma_res_to_ttop_list function
#' 
#'
#' @param ttop_list  list of limma top table results obtained from limma_res_to_ttop_list
#' @return data.frame of t-values, where columns are the contrasts and rows are the omic features tested
#' @export
ttop_list_to_t_table <- function(ttop_list)
{
  if(length(ttop_list) > 1)
  {
    t_table <- merge(ttop_list[[1]][,c(1,4)], ttop_list[[2]][,c(1,4)], by = "ID")
    if(length(ttop_list) > 2)
    {
      for(i in 3:length(ttop_list))
      {
        t_table <- merge(t_table, ttop_list[[i]][,c(1,4)], by = "ID")
      }
    }
  }
  else
  {
    t_table <- ttop_list[[1]][,c(1,4)]
  }
  names(t_table) <- c("ID",names(ttop_list))
  return(t_table)
}

#' ttop_list_to_log2FC_table
#'
#' same as ttop_list_to_t_table, but with log2FC instead of t-values
#'
#' @param ttop_list  list of limma top table results obtained from limma_res_to_ttop_list
#' @return data.frame of t-values, where columns are the contrasts and rows are the omic features tested
#' @export
ttop_list_to_log2FC_table <- function(ttop_list)
{
  if(length(ttop_list) > 1)
  {
    log2FC_table <- merge(ttop_list[[1]][,c(1,2)], ttop_list[[2]][,c(1,2)], by = "ID")
    if(length(ttop_list) > 2)
    {
      for(i in 3:length(ttop_list))
      {
        log2FC_table <- merge(log2FC_table, ttop_list[[i]][,c(1,2)], by = "ID")
      }
    }
  }
  else
  {
    log2FC_table <- ttop_list[[1]][,c(1,2)]
  }
  names(log2FC_table) <- c("ID",names(ttop_list))
  return(log2FC_table)
}

#' limma_res_to_ttop_list
#'
#' This function converts the output of the runLimma function into a list of
#' limma top tables
#'
#' @param limma_res  output of the runLimma function
#' @param comp_names  identifiers of the contrasts that were considered in the 
#' runLimma function
#' @param number  number of top features to be included in the top table results
#' @param adjust.method  method to adjust p-values, same as for the topTable limma function
#' @return a list of limma top table results
#' @export
limma_res_to_ttop_list <- function(limma_res, comp_names, number, adjust.method = "fdr")
{
  ttop_list <- list()
  n_comp <- length(limma_res[[2]][1,])
  for(i in 1:n_comp)
  {
    ttop_list[[i]] <- ttopFormatter(topTable(limma_res[[1]], coef = i, number = number, adjust.method = adjust.method))
    ttop_list[[i]] <- ttop_list[[i]][complete.cases(ttop_list[[i]]),]
  }
  names(ttop_list) <- comp_names
  return(ttop_list)
}

#' t_table_metactivity_input_formater
#'
#' This function allows to format a table of t-values (or log2FCs) for the 
#' metactivity function. Essentially, it adds a compartment code (using BIGG
#' nomenclature) at the end of each metabolite identifier, so that they can be 
#' mapped on the compartmentalized reaction prior knowledge network.
#'
#' @param metabolomic_t_table  the metabolomic t-table (or log2FC table) that 
#' is returned by ttop_list_to_t_table (or ttop_list_to_log2FC_table)
#' @param mapping_table  a mapping table to translate the metabolite 
#' identifiers of the t-table into KEGG IDs
#' @param affixes  a character vectors corresponding to the compartment affixes
#' to be considered. If you wish to only focus on one compartment, 
#' e.g. mitochondria, you can input only c("m") for example.
#' @return a formatted t-table (or log2FC table) for use with the metactivity function
#' @export
t_table_metactivity_input_formater <- function(metabolomic_t_table, mapping_table, affixes = c("c","l","x","m","e","n","r"))
{
  names(metabolomic_t_table)[1] <- "metabolite"
  metabolomic_t_table[,1] <- gsub(",","_",metabolomic_t_table[,1])
  metabolomic_t_table[,1] <- gsub(" ","",metabolomic_t_table[,1])

  names(mapping_table)[1] <- "metabolite"
  mapping_table$metabolite <- gsub(",","_",mapping_table$metabolite)
  mapping_table$metabolite <- gsub(" ","",mapping_table$metabolite)

  if(!grepl("cpd:",mapping_table[1,"KEGG"]))
  {
    mapping_table$KEGG <- paste("cpd:",mapping_table$KEGG, sep = "")
  }

  temp_mapping_table <- mapping_table

  for (aff in affixes)
  {
    new_mapping_table <- temp_mapping_table
    new_mapping_table$KEGG <- paste(new_mapping_table$KEGG, aff, sep = "_")
    mapping_table <- as.data.frame(rbind(mapping_table, new_mapping_table))
  }

  mapping_table[,1] <- tolower(mapping_table[,1])
  metabolomic_t_table$metabolite <- tolower(metabolomic_t_table$metabolite)

  metabolomic_t_table <- merge(metabolomic_t_table, mapping_table, by = "metabolite")
  table_length <- length(metabolomic_t_table[1,])
  metabolomic_t_table <- metabolomic_t_table[,c(table_length,2:(table_length-1))]

  return(metabolomic_t_table)
}

#'\code{ttopFormatter}
#'
#'This function is simply designed to format the toptable of limma with first column as gene identifiers instead of only row.names.
#'
#'@param ttop a toptable dataframe generated by the topTable function of limma
#'
#'@return a dataframe similar to the output of topTable function of limma but with first column as IDs instead of onyl row.names.
ttopFormatter <- function(ttop)
{
  ttop$ID <- row.names(ttop)
  ttop <- ttop[,c(7,1,2,3,4,5,6)]
  ttop <- ttop[complete.cases(ttop),]
  return(ttop)
}

#'\code{nicePCA}
#'
#'This function generate a 3*3 arrangeGrob plot object (that can be subsequently diplayed or saved).
#'Each cell of the 3*3 plot grid correspond to a specific representation of the result of a principal component analysis performed on a measurment dataframe.
#'The first input is a n*m data.frame, where n is the number of measured omic features (genes, proteins, metabolites...) and m is the number of samples.
#'The second input is a basic n*2 target dataframe (such as generated by the generateTarget function), where n is the number of samples.
#'
#' @param df the measurment n*m dataframe (n is number of omic features, m is number of samples) where columns are ordered by conditions.
#' @param targets A n*2 dataframe, where n is the number of samples. First column correspond to samples, second column correspond to conditions.
#' @param components a vector of three integers, corresponding to the components to be plotted
#' @param centering a boolean parameter to indicate wether samples should be mean centered
#' @param scaling a boolean parameter to indicate wether samples should be scaled (x/variance)
#' @param pointSize an integer parameter to indicate the desired point size for the components scatter plots.
#' @return an 3*3 arrangeGrob object containing various graphical representation of the result of a PCA.
#' @import ggplot2
#' @import grid
#' @importFrom gridExtra arrangeGrob
#' @importFrom cowplot get_legend
#' @export
nicePCA <- function(df, targets, components = c(1,2,3), centering  = T, scaling = F, pointSize = 4, no_label = FALSE)
{  
  if (!is.null(targets))
  {
    df_and_targets <- make_df_and_targets_great_again(df,targets)
    df <- df_and_targets[[1]]
    targets <- df_and_targets[[2]]
  }
  
  sample <- names(df)
  pca <- prcomp(t(df[complete.cases(df),]), center = centering, scale. = scaling)
  explained <- (pca$sdev)^2 / sum(pca$sdev^2)
  
  xCompLab <- paste(paste("PC",components[1], sep = "")," (", sep = "")
  yCompLab <- paste(paste("PC",components[2], sep = "")," (", sep = "")
  zCompLab <- paste(paste("PC",components[3], sep = "")," (", sep = "")
  
  data.to.plot <- as.data.frame(cbind(pca$x[,components[1]], pca$x[,components[2]], pca$x[,components[3]]))
  if ("color" %in% names(targets))
  {
    data.to.plot$condition <- targets$color
  }
  else
  {
    data.to.plot$condition <- targets$condition
  }
  data.to.plot$sample <- targets$sample
  colnames(data.to.plot) <- c("pc1", "pc2", "pc3","condition", "sample")
  
  comp1 <- ggplot(data.to.plot, aes(x = pc1)) + geom_density(aes(alpha = 0.3)) + theme_minimal() +
    xlab(paste(xCompLab,round(100*explained[components[1]],digits=2),'%)',sep='')) +
    theme(legend.position = "none") +
    scale_x_continuous(position = "top", limits = c(-max(abs(pca$x[,components[1]])),max(abs(pca$x[,components[1]])))) +
    scale_y_continuous(position = "left")
  comp2 <- ggplot(data.to.plot, aes(x = pc2)) + geom_density(aes(alpha = 0.3)) + theme_minimal() +
    xlab(paste(yCompLab,round(100*explained[components[2]],digits=2),'%)',sep='')) +
    theme(legend.position = "none") +
    scale_x_continuous(position = "top", limits = c(-max(abs(pca$x[,components[2]])),max(abs(pca$x[,components[2]])))) +
    scale_y_continuous(position = "right")
  comp3 <- ggplot(data.to.plot, aes(x = pc3)) + geom_density(aes(alpha = 0.3)) + theme_minimal() +
    xlab(paste(zCompLab,round(100*explained[components[3]],digits=2),'%)',sep='')) +
    theme(legend.position = "none") +
    scale_x_continuous(position = "top", limits = c(-max(abs(pca$x[,components[3]])),max(abs(pca$x[,components[3]])))) +
    scale_y_continuous(position = "right")
  
  just_for_legend <- ggplot(data.to.plot, aes(x=pc1, y=pc2, color=condition)) +
    geom_point(size=pointSize, alpha = 0.3) +
    geom_text(aes(label=sample), nudge_y = 1.0, size=6) +
    scale_alpha_discrete(range=c(0.3, 1.0)) +
    #geom_path(arrow=arrow()) +
    theme_minimal() +
    xlab(paste(xCompLab,round(100*explained[components[1]],digits=2),'%)',sep='')) +
    ylab(paste(yCompLab,round(100*explained[components[2]],digits=2),'%)',sep='')) +
    xlim(c(-max(abs(pca$x[,components[1]])),max(abs(pca$x[,components[1]])))) +
    ylim(c(-max(abs(pca$x[,components[2]])),max(abs(pca$x[,components[2]]))))
  
  legend <- get_legend(just_for_legend)
  rm(just_for_legend)
  
  if (no_label)
  {
    comp1vcomp2 <- ggplot(data.to.plot, aes(x=pc1, y=pc2, color=condition)) +
      geom_point(size=pointSize, alpha = 0.3) +
      scale_alpha_discrete(range=c(0.3, 1.0)) +
      #geom_path(arrow=arrow()) +
      theme_minimal() +
      xlab("") +
      ylab(paste(yCompLab,round(100*explained[components[2]],digits=2),'%)',sep='')) +
      ylim(c(-max(abs(pca$x[,components[2]])),max(abs(pca$x[,components[2]])))) + theme(legend.position = "none") +
      scale_x_continuous(position = "top",limits = c(-max(abs(pca$x[,components[1]])),max(abs(pca$x[,components[1]]))))
    
    comp1vcomp3 <- ggplot(data.to.plot, aes(x=pc1, y=pc3, color=condition)) +
      geom_point(size=pointSize, alpha = 0.3) +
      scale_alpha_discrete(range=c(0.3, 1.0)) +
      #geom_path(arrow=arrow()) +
      theme_minimal() +
      xlab("") +
      ylab(paste(zCompLab,round(100*explained[components[3]],digits=2),'%)',sep='')) +
      ylim(c(-max(abs(pca$x[,components[3]])),max(abs(pca$x[,components[3]])))) + theme(legend.position = "none") +
      scale_x_continuous(position = "top", limits = c(-max(abs(pca$x[,components[1]])),max(abs(pca$x[,components[1]]))))
    
    comp2vcomp3 <- ggplot(data.to.plot, aes(x=pc2, y=pc3, color=condition)) +
      geom_point(size=pointSize, alpha = 0.3) +
      scale_alpha_discrete(range=c(0.3, 1.0)) +
      #geom_path(arrow=arrow()) +
      theme_minimal() +
      xlab("") +
      ylab("") +
      scale_y_continuous(position = "right", limits = c(-max(abs(pca$x[,components[3]])),max(abs(pca$x[,components[3]])))) +
      scale_x_continuous(position = "top", limits = c(-max(abs(pca$x[,components[2]])),max(abs(pca$x[,components[2]])))) + theme(legend.position = "none")
  }
  else
  {
    comp1vcomp2 <- ggplot(data.to.plot, aes(x=pc1, y=pc2, color=condition)) +
      geom_point(size=pointSize, alpha = 0.3) +
      geom_text_repel(aes(label=sample), size = 2) +
      scale_alpha_discrete(range=c(0.3, 1.0)) +
      #geom_path(arrow=arrow()) +
      theme_minimal() +
      xlab("") +
      ylab(paste(yCompLab,round(100*explained[components[2]],digits=2),'%)',sep='')) +
      ylim(c(-max(abs(pca$x[,components[2]])),max(abs(pca$x[,components[2]])))) + theme(legend.position = "none") +
      scale_x_continuous(position = "top",limits = c(-max(abs(pca$x[,components[1]])),max(abs(pca$x[,components[1]]))))
    
    comp1vcomp3 <- ggplot(data.to.plot, aes(x=pc1, y=pc3, color=condition)) +
      geom_point(size=pointSize, alpha = 0.3) +
      geom_text_repel(aes(label=sample), size = 2) +
      scale_alpha_discrete(range=c(0.3, 1.0)) +
      #geom_path(arrow=arrow()) +
      theme_minimal() +
      xlab("") +
      ylab(paste(zCompLab,round(100*explained[components[3]],digits=2),'%)',sep='')) +
      ylim(c(-max(abs(pca$x[,components[3]])),max(abs(pca$x[,components[3]])))) + theme(legend.position = "none") +
      scale_x_continuous(position = "top", limits = c(-max(abs(pca$x[,components[1]])),max(abs(pca$x[,components[1]]))))
    
    comp2vcomp3 <- ggplot(data.to.plot, aes(x=pc2, y=pc3, color=condition)) +
      geom_point(size=pointSize, alpha = 0.3) +
      geom_text_repel(aes(label=sample), size = 2) +
      scale_alpha_discrete(range=c(0.3, 1.0)) +
      #geom_path(arrow=arrow()) +
      theme_minimal() +
      xlab("") +
      ylab("") +
      scale_y_continuous(position = "right", limits = c(-max(abs(pca$x[,components[3]])),max(abs(pca$x[,components[3]])))) +
      scale_x_continuous(position = "top", limits = c(-max(abs(pca$x[,components[2]])),max(abs(pca$x[,components[2]])))) + theme(legend.position = "none")
    
  }
  
  #barplot(explained, col = "lightblue")
  explained <- as.data.frame(explained)
  explained$Components <- as.numeric(row.names(explained))
  names(explained)[1] <- "% of explained variance"
  
  boulder <- ggplot(explained, aes(x = Components, y = `% of explained variance` , fill = "red")) +
    geom_col() +
    theme_minimal() +
    scale_y_continuous(position = "right") +
    theme(legend.position = "none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_x_continuous(position = "top", breaks = c(1:length(explained[,1])))
  
  return(arrangeGrob(comp1, legend, boulder,comp1vcomp2, comp2, legend, comp1vcomp3, comp2vcomp3, comp3, ncol = 3))
}

#'\code{make_df_and_targets_great_again}
#'
#'This function Check wether sample names are coherent between the measurment dataframe and the target dataframe. If they are coherent,
#' the target dataframe rows are reordered to match the column order of the measurment dataframe.
#'@param df the measurment n*m dataframe (n is number of omic features, m is number of samples) where columns are ordered by conditions.
#'@param targets A n*2 dataframe, where n is the number of samples. First column correspond to samples, second column correspond to conditions.
#'@return A list, first element is the cleaned measurment dataframe and second element is the cleaned target dataframe.
make_df_and_targets_great_again <- function(df, targets)
{
  names(targets)[c(1,2)] <- c("sample","condition")
  bad_samples <- c(NA)
  i <- 1
  for (sample in names(df))
  {
    if (!(sample %in% targets[,1]))
    {
      bad_samples <- c(bad_samples, sample)
    }
    else
    {
      if (i == 1)
      {
        clean_targets <- as.data.frame(matrix(NA, 1, length(targets[1,])))
        names(clean_targets) <- names(targets)
        clean_targets[1,] <- targets[targets$sample == sample,]
        i <- i+1
      }
      else
      {
        clean_targets <- as.data.frame(rbind(clean_targets,targets[targets$sample == sample,]))
        i <- i+1
      }
    }
  }
  targets <- clean_targets
  if (length(bad_samples) > 1)
  {
    bad_samples <- bad_samples[-1]
    print(paste("These samples were not found : ", bad_samples, sep = ""))
    print(bad_samples)
    df <- df[,!(names(df) %in% bad_samples)]
  }
  return(list(df,targets))
}
