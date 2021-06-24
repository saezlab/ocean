# These are generated from the tutorial, but didn't want to make any changes there, 
# so that's why I save them here
# Using the tutorial script to save the expected (test) outputs would be much better :)

ocean_path <- system.file(package="ocean")

# test_limma
unique(toy_targets$condition)
comparisons <- list('tumorVsHealthy' = c(1,-2))
limmaRes <- runLimma(measurements = toy_metabolomic_data,
                     targets = toy_targets,
                     comparisons = comparisons)
saveRDS(limmaRes, file.path(ocean_path, "testdata/limma_res.rds"))

# test_ttop
t_table <- ttop_list_to_t_table(
  limma_res_to_ttop_list(limma_res = limmaRes,
                         comp_names = names(comparisons),
                         number = length(toy_metabolomic_data[,1]),
                         adjust.method = "fdr"))
saveRDS(t_table, file.path(ocean_path, "testdata/ttop_res.rds"))


# test_metabolic_mapp
t_table <- t_table_metactivity_input_formater(metabolomic_t_table = t_table,
                                              mapping_table = mapping_table,
                                              affixes = c("c","l","x","m","e","n","r"))
saveRDS(t_table, file.path(ocean_path, "testdata/meta_map.rds"))



penalty_min <- 6 #minimum 1 and integer
penalty_max <- 6 #maximum 9 and integer

# Test network build and filter
TCA_network_nocofact <-
  model_to_pathway_sif(pathway_to_keep = c("Citric acid cycle",
                                           "Glycolysis/gluconeogenesis",
                                           "Pyruvate metabolism",
                                           "Transport, mitochondrial")) %>%
  remove_cofactors() %>%
  compress_transporters() %>%
  split_transaminases()
saveRDS(TCA_network_nocofact, file.path(ocean_path, "testdata/network_filt.rds"))



# Test Forest
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
saveRDS(regulons_df, file.path(ocean_path, "testdata/met_regulons.rds"))


# Test metabolic activity scores
metactivity_res <- metactivity(metabolomic_t_table = t_table,
                               regulons_df = regulons_df, 
                               compartment_pattern = "_[a-z]$", 
                               k = 1000)

saveRDS(metactivity_res, file.path(ocean_path, "testdata/met_activity.rds"))

mean_ES_df <- metactivity_res$ES
mean_NES_df <- metactivity_res$NES

translated_results <- translate_results(regulons_df = regulons_df,
                                        t_table = t_table,
                                        mapping_table = mapping_table)
saveRDS(translated_results, file.path(ocean_path, "testdata/translated_res.rds"))


# Test Enzyme Plots
plots <- plotMetaboliteContribution(enzyme = '4967_1738_8050_1743', 
                                    stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_results$regulons_df, 
                                    contrast_index = 1, 
                                    stat_name = 't', 
                                    scaling_factor = 5, nLabels = 15)
saveRDS(plots, file.path(ocean_path, "testdata/enzyme_plots.rds"))


# Test Pathway Plot
hm <- pathway_HM(mean_NES_df = mean_NES_df,
                 pathway_name = 'KEGG_CITRATE_CYCLE_TCA_CYCLE',
                 pathways = kegg_pathways)
saveRDS(hm, file.path(ocean_path, "testdata/pathway_plot.rds"))
