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

#' makeContrastsAlt
#'
#' ipsum...
#'
#' @param targets ipsum...
#' @param comparisons ipsum...
#' @return ipsum...
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

#' runLimma
#'
#' ipsum...
#'
#' @param measurments  ipsum...
#' @param targets ipsum...
#' @param comparisons ipsum...
#' @param regress_out ipsum...
#' @return ipsum...
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
#' ipsum...
#'
#' @param ttop_list  ipsum...
#' @return ipsum...
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
#' ipsum...
#'
#' @param ttop_list  ipsum...
#' @return ipsum...
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
#' ipsum...
#'
#' @param limma_res  ipsum...
#' @param comp_names  ipsum...
#' @param number  ipsum...
#' @param adjust.method  ipsum...
#' @return ipsum...
#' @export
limma_res_to_ttop_list <- function(limma_res, comp_names, number, adjust.method = "fdr")
{
  ttop_list <- list()
  n_comp <- length(limma_res[[2]][1,])
  for(i in 1:n_comp)
  {
    ttop_list[[i]] <- ttopFormatter(topTable(limmaRes[[1]], coef = i, number = number, adjust.method = adjust.method))
    ttop_list[[i]] <- ttop_list[[i]][complete.cases(ttop_list[[i]]),]
  }
  names(ttop_list) <- comp_names
  return(ttop_list)
}

#' t_table_metactivity_input_formater
#'
#' ipsum...
#'
#' @param metabolomic_t_table  ipsum...
#' @param mapping_table  ipsum...
#' @param affixes  ipsum...
#' @return ipsum...
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

#' ttopFormatter
#'
#' ipsum...
#'
#' @param ttop  ipsum...
#' @return ipsum...
ttopFormatter <- function(ttop)
{
  ttop$ID <- row.names(ttop)
  ttop <- ttop[,c(7,1,2,3,4,5,6)]
  ttop <- ttop[complete.cases(ttop),]
  return(ttop)
}
