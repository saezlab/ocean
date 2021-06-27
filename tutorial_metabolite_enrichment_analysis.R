################ ocEAn Tutorial Metabolite Enrichment Analysis #################

##Install ocEAn
library(devtools)
install_github("saezlab/ocean")


##Tutorial with a kidney cancer toy metabolomic dataset
library(ocean)
library(decoupleR)


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

##Select all pathways from recon2_redhuman$pathway to be included in network
all_pathways <- as.vector(unique(recon2_redhuman$pathway))
all_pathways <- all_pathways[!is.na(all_pathways)] 

##Create a reaction network from the recon_redhuman model (metabolites & enzymes)
reaction_network <- model_to_pathway_sif(pathway_to_keep = all_pathways)

##Translate enzyme complexes by mapping identifiers to names
reaction_network <- translate_complexes(reaction_network)

##Create data frame metabolite - affiliated enzyme instead of source-target
metabolites <- rearrange_dataframe(reaction_network)

##Remove compartment information and "cpd:" to get pure KEGG IDs
metabolites$metabolites <- get_pure_kegg_ids(metabolites$metabolites)
metabolites <- distinct(metabolites)  #keep only unique rows

##Map pathways to metabolites
metabolites_pathway_df <- map_pathways_to_metabolites(metabolites)

#Save metabolites-pathways data frame
#write.csv(metabolites_pathway_df,"./results/metabolites_pathway_df.csv")


##Metabolite enrichment analysis with package decoupleR

#'Calculate the activity enrichment score and p-value of all pathways by using
#'the conditions in the matrix (patient samples vs metabolite expression) by
#'calculating the mean over the expression of all metabolites.

##Prepare input data
network <- metabolites_pathway_df    #input 1: pathways and their associated metabolites
network$mor <- 1                     #add new column for mode of regulation
network$likelihood <- 1              #add new column for edge likelihood

t_table_kegg <- t_table
t_table_kegg$KEGG <- get_pure_kegg_ids(t_table_kegg$KEGG) #remove compartment info
t_table_kegg <- unique(t_table_kegg)
row.names(t_table_kegg) <- t_table_kegg$KEGG  #convert kegg column to row names
t_table_kegg$KEGG <- NULL            #delete kegg ids column
t_table_kegg <- unique(t_table_kegg) #ensure only 1 kegg is mapped to 1 metabolite
mat <- as.matrix(t_table_kegg)       #input 2: expression matrix

##Perform enrichment analysis
enrichment <- run_mean(mat, network,
                          .source = .data$pathway,
                          .target = .data$metabolites,
                          .mor = .data$mor, 
                          .likelihood = .data$likelihood,  
                       times = 10000,            #number of permutations
                       seed = 42,                #a single integer
                       sparse = TRUE,
                       randomize_type = "rows")  #randomize matrix
enrichment$condition <- NULL
enrichment <- as.data.frame(enrichment)
colnames(enrichment) <- c("statistic", "pathway", "score", "p-value")

##Keep only rows with statistic "normalized_mean"
enrichment_norm <- enrichment[enrichment$statistic != "mean", ]

##Visualize most significant pathways in a bar plot
score <- 1.0  #define cutoff score
barplot = plot_significant_pathways(enrichment_norm, score)
barplot

#ggsave("enriched_pathways.png", plot = barplot,
#       path = "results/", scale = 1, dpi = 300, limitsize = TRUE)
