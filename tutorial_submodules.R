library(ocean)

##Differential analysis
unique(toy_targets$condition)
comparisons <- list('tumorVsHealthy' = c(1,-2))

limmaRes <- runLimma(measurements = toy_metabolomic_data, targets = toy_targets, comparisons = comparisons)

##Format differeential anaylsis result
t_table <- ttop_list_to_t_table(
  limma_res_to_ttop_list(limma_res = limmaRes,
                         comp_names = names(comparisons),
                         number = length(toy_metabolomic_data[,1]),
                         adjust.method = "fdr"))

##This step is particularly important because this is where the users metabolic identifiers are mapped to the kegg ids used by the method.
#Thus, the user should provide a mapping table in the same format as the mapping_table presented here (you can look at it to inspire yourself from it)
#The mapping table should map the users own metabolic identifiers to kegg compound IDs
t_table <- t_table_metactivity_input_formater(metabolomic_t_table = t_table,
                                              mapping_table = mapping_table,
                                              affixes = c("c","l","x","m","e","n","r"))

##Prepare the metabolic enzyme sets
penalty_min <- 6 #minimum 1 and integer
penalty_max <- 6 #maximum 9 and integer

## These are the available pathways to choose from
View(unique(recon2_redhuman$pathway))

##Select pathways relevant to include in network
TCA_network <- model_to_pathway_sif(pathway_to_keep = c("Citric acid cycle",
                                                        "Glycolysis/gluconeogenesis",
                                                        "Pyruvate metabolism",
                                                        "Transport, mitochondrial"))

TCA_network_nocofact <- remove_cofactors(TCA_network)

##This is to simplify the network sutrcture by compressing redundant transporters
TCA_network_nocofact <- compress_transporters(sub_network_nocofact = TCA_network_nocofact)

##This is to avoid cross ping pong between reactants and products of reversible transaminases 
TCA_network_nocofact <- split_transaminases(sub_network_nocofact = TCA_network_nocofact)

enzymes <- unique(TCA_network_nocofact$attributes$V1)
enzymes <- enzymes[!grepl("_[clxmenr]$",enzymes)]

#branch_length applies a cutoff on the minimum length of the reaction network 
#upstream and down stream of a given enzyme. Here we use a minimum length of 3 for
#both directions
TCA_forest <- forestMaker(enzymes, TCA_network_nocofact$reaction_network, branch_length = c(3,3))

reaction_set_list <- prepare_metabolite_set(penalty_range = penalty_min:penalty_max,  
                                            forest = TCA_forest,
                                            measured_metabolites = t_table$KEGG)

reaction_set_list_merged <- condense_metabolite_set(reaction_set_list = reaction_set_list)

penalty <- 6 #has be between penalty_min and penalty_max and integer

regulons_df <- prepare_regulon_df(reaction_set_list_merged, penalty, c(0,1))

##Compute metabolic enzme enrichment score
metactivity_res <- metactivity(metabolomic_t_table = t_table, 
                               regulons_df = regulons_df, 
                               compartment_pattern = "_[a-z]$", 
                               k = 1000)

mean_ES_df <- metactivity_res$ES

mean_NES_df <- metactivity_res$NES

##translate the metabolic ids back to names
translated_results <- translate_results(regulons_df = regulons_df, t_table = t_table, mapping_table = mapping_table)

##Visualise results for single enzmes
plots <- plotMetaboliteContribution(enzyme = '4967_1738_8050_1743', 
                                    stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_results$regulons_df, 
                                    contrast_index = 1, 
                                    stat_name = 't', 
                                    scaling_factor = 5, nLabels = 15)

plot(plots$scatter)
plot(plots$cumsumPlot)

##Visualise results at KEGG pathway level
hm <- pathway_HM(mean_NES_df = mean_NES_df, pathway_name = 'KEGG_CITRATE_CYCLE_TCA_CYCLE', pathways = kegg_pathways)
plot(hm)

## visualise the network
plot_reaction_network(TCA_network_nocofact, t_table, mean_NES_df, column_index = 1, vis.height = 2000)

