########################## ocEAn Tutorial ORA ##################################

## install ocEAn
#library(devtools)
#install_github("saezlab/ocean")


## Tutorial with a kidney cancer toy metabolomic dataset
library(ocean)
library(piano)
source("./R/translate_complexes.R")
source("./R/support_metabolite_ORA.R")


##Differential analysis
unique(toy_targets$condition)
comparisons <- list('tumorVsHealthy' = c(1,-2))

limmaRes <- runLimma(measurements = toy_metabolomic_data,
                     targets = toy_targets,
                     comparisons = comparisons)

##Format differential analysis result
t_table <- ttop_list_to_t_table(
  limma_res_to_ttop_list(limma_res = limmaRes,
                         comp_names = names(comparisons),
                         number = length(toy_metabolomic_data[,1]),
                         adjust.method = "fdr"))

t_table <- t_table_metactivity_input_formater(metabolomic_t_table = t_table,
                                              mapping_table = mapping_table,
                                              affixes = c("c","l","x","m","e","n","r"))


## These are the available pathways to choose from
View(unique(recon2_redhuman$pathway))

##Select pathways relevant to include in network
TCA_network <- model_to_pathway_sif(pathway_to_keep = c("Citric acid cycle",
                                                        "Glycolysis/gluconeogenesis",
                                                        "Pyruvate metabolism",
                                                        "Transport, mitochondrial"))

##Translate enzyme complexes by mapping identifiers to names
TCA_network_trans <- translate_complexes(TCA_network)


##Create data frame metabolite - affiliated enzyme instead of source-target
metabolites <- rearrange_dataframe(TCA_network_trans)

#rename metabolites column for translation
colnames(metabolites)[colnames(metabolites) == "metabolites"] <- "targets" 

##Translate the metabolic ids back to names
translated_results <- translate_results(regulons_df = metabolites,
                                        t_table = t_table,
                                        mapping_table = mapping_table)

names(translated_results) <- c("t_table", "metabolites_trans")
colnames(translated_results$t_table)[colnames(translated_results$t_table) == "KEGG"] <- "metabolites"  
colnames(translated_results$metabolites_trans)[colnames(translated_results$metabolites_trans) == "targets"] <- "metabolites"


##Remove compartment information
for (i in 1:length(translated_results)) {

  translated_results[[i]]["metabolites"] <- sapply(translated_results[[i]]["metabolites"],
                                                   function(x){
                                                    x <- gsub("_.$","",x)
                                                    return(x)
                                                   }, simplify = F, USE.NAMES = F)
  translated_results[[i]] <- distinct(translated_results[[i]])
}


##Map pathways to metabolites
metabolites_pathway_df <- map_pathways_to_metabolites(translated_results$metabolites_trans)

# save output files as RData
saveRDS(metabolites_pathway_df,"./results/metabolites_pathway_df.RData")


##Over Representation Analysis (ORA) with metabolites 

#'The ORA is performed to determine whether metabolites with statistically
#'significant differential abundances interact preferentially with metabolic
#'enzymes belonging to the same metabolic pathways of interest. The interaction
#'of significant metabolites with pathways is statistically tested by comparing
#'metabolites with differential abundances to all input metabolites.
#'
#'The results of the ORA might help to elucidate significant differences between
#'two biological conditions such as cancer vs. healthy tissue, or treatment vs.
#'no treatment. The metabolites with unusual abundances give information about
#'which metabolic pathways are SIGNIFICANTLY different from what is expected.
#'This information might be helpful to identify metabolic enzymes that are
#'presumably deregulated.


##Create a vector of all input metabolites (background)
#translated metabolites from t-table are used because names from input could be slightly different
metabolites_universe <- translated_results$t_table$metabolites 
#metabolites_universe <- unlist(unique(metabolites_pathway_df$metabolites))  # all input metabolites belonging to pathways


##Create a vector of all significant metabolites
p_value <- 2.0 #define cutoff p-value
metabolites_significant_df <- translated_results$t_table[translated_results$t_table$tumorVsHealthy <= -p_value | translated_results$t_table$tumorVsHealthy >= p_value, ] 
metabolites_significant <- unique(metabolites_significant_df$metabolites)


##We perform the ORA with Fisher's exact test (from piano package)
sig_pathways <- runGSAhyper(genes = metabolites_significant, 
                            universe = metabolites_universe,
                            gsc = loadGSC(metabolites_pathway_df),
                            adjMethod = "fdr")
sig_pathways_df <- as.data.frame(sig_pathways$resTab)  %>% 
  tibble::rownames_to_column(var = "pathway")



##Visualize results


