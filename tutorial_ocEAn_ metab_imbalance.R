################################# OCEAN ########################################

#'R package for metabolic enzyme enrichment analysis
#'
#'OCEAN: a method that defines metabolic enzyme footprint from a curated reduced
#'       version of the recon2 reaction network. The metabolic enzyme footprints
#'       are used to explore coordinated deregulations of metabolite abundances
#'       with respect to their position relative to metabolic enzymes.
#'       This is similar to Kinase-substrate and TF-targets enrichment analyses.

## install ocEAn
library(devtools)
install_github("saezlab/ocean")


## Tutorial with a kidney cancer toy metabolomic dataset
library(ocean)

##Differential analysis 
#with data from: https://www.embopress.org/doi/full/10.15252/msb.20209730

unique(toy_targets$condition)
comparisons <- list('tumorVsHealthy' = c(1,-2))

limmaRes <- runLimma(measurements = toy_metabolomic_data,
                     targets = toy_targets,
                     comparisons = comparisons)

#this differential anaylsis result represents metabolic up and down-regulations
#of metabolite abundances in kidney tumor tissue compared to adjacent 
#non-tumoral tissue

##Format differential analysis result

#In order to assess the metabolic imbalance, we use the t statistic of the 
#differential analysis. The t-statistic represents a good estimate of the 
#relative significance of up and down-regulations of metabolic abundances.

t_table <- ttop_list_to_t_table(
  limma_res_to_ttop_list(limma_res = limmaRes,
                         comp_names = names(comparisons),
                         number = length(toy_metabolomic_data[,1]),
                         adjust.method = "fdr"))

##Next, the metabolic names have to be converted to a standard identifer.
#In this case, we use KEGG IDs. 
#This step is particularly important because this is where the users 
#metabolic identifiers are mapped to the kegg ids used by the method.

#THE USER SHOULD PROVED A MAPPING TABLE IN THE SAME FORMAT AS THE ONE PRESENTED HERE
#(you can look at it to inspire yourself from it).
#The mapping table should map the users own metabolic identifiers to kegg compound IDs
t_table <- t_table_metactivity_input_formater(metabolomic_t_table = t_table,
                                              mapping_table = mapping_table,
                                              affixes = c("c","l","x","m","e","n","r"))


##These are the available pathways to choose from
View(unique(recon2_redhuman$pathway))

######## SUBNETWORKS
data(expressed_genes)

#In this case, we will use all pathways available. 
#Users can restrict this list on a per-need basis.

all_pathways <- unique(recon2_redhuman$pathway)
sub_network <- model_to_pathway_sif(pathway_to_keep = all_pathways$X1)

sub_network <- translate_complexes(sub_network)

sub_network_nocofact <- remove_cofactors(sub_network)

#Filter out any gene that isn't expressed from the metabolic network
non_expressed_genes <- expressed_genes[rowSums(expressed_genes[,c(2),drop = F]) == 0,c(1),drop = F] 
non_expressed_genes <- c(non_expressed_genes$gene)

tokeep <- list()
for(i in 1:length(sub_network_nocofact$reaction_network[,1]))
{
  elements <- gsub("[><].*","",sub_network_nocofact$reaction_network[i,])
  elements <- elements[!grepl("cpd:",elements)]
  elements <- unlist(strsplit(elements, "_"))
  tokeep[[i]] <- sum(elements %in% non_expressed_genes) == 0
}
tokeep <- unlist(tokeep)
sub_network_nocofact$reaction_network <- sub_network_nocofact$reaction_network[tokeep,]

#Compress transporters to simplify the network
sub_network_nocofact <- compress_transporters(sub_network_nocofact = sub_network_nocofact)

#Split transaminases to preserve correct metabolic transformation routes
sub_network_nocofact <- split_transaminases(sub_network_nocofact = sub_network_nocofact)

enzymes <- unique(sub_network_nocofact$attributes$V1)
enzymes <- enzymes[!grepl("_[clxmenr]$",enzymes)]

#Convert the network into a forest (list of enzymes (trees) with correspoding metabolic signatures (branches))
sub_forest <- forestMaker(enzymes, sub_network_nocofact$reaction_network, branch_length = c(1,1), remove_reverse = T)

###################
##Prepare the metabolic enzyme sets
penalty_min <- 6 #minimum 1 and integer
penalty_max <- 8 #maximum 9 and integer

reaction_set_list <- prepare_metabolite_set(penalty_range = penalty_min:penalty_max,   
                                            forest = sub_forest,
                                            measured_metabolites = t_table$KEGG)

reaction_set_list_merged <- condense_metabolite_set(reaction_set_list = reaction_set_list)

penalty <- 8 #has be between penalty_min and penalty_max and integer

regulons_df <- prepare_regulon_df(reaction_set_list_merged, penalty, filter_imbalance = c(0,1))

##Compute metabolic enzme enrichment score
metactivity_res <- metactivity(metabolomic_t_table = t_table, 
                               regulons_df = regulons_df, 
                               compartment_pattern = "_[a-z]$", 
                               k = 1000)

mean_ES_df <- metactivity_res$ES

mean_NES_df <- metactivity_res$NES

##########################
translated_results <- translate_results(regulons_df = regulons_df, t_table = t_table, mapping_table = mapping_table, compress_compartments = T)

##Visualise results for single enzmes
plots <- plotMetaboliteContribution(enzyme = 'OGDH_DLD_PDH_DLST>1052', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_results$regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 6, nLabels =  25)

plot(plots$scatter)

plots <- plotMetaboliteContribution(enzyme = 'BCAT1', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_results$regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 6, nLabels =  25)

plot(plots$scatter)
