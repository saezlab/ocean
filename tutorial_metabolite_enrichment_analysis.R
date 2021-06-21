################ ocEAn Tutorial Metabolite Enrichment Analysis #################

## install ocEAn
library(devtools)
install_github("saezlab/ocean")


## Tutorial with a kidney cancer toy metabolomic dataset
library(ocean)
library(decoupleR)
library(piano)


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

##Create network with metabolites/enzymes as sources/targets
network <- model_to_pathway_sif(pathway_to_keep = all_pathways)

##Translate enzyme complexes by mapping identifiers to names
network_trans <- translate_complexes(network)

##Create data frame metabolite - affiliated enzyme instead of source-target
metabolites <- rearrange_dataframe(network_trans)

##Remove compartment information and "cpd:" to get pure KEGG IDs
metabolites$metabolites <- get_pure_kegg_ids(metabolites$metabolites)

##Keep only unique rows
metabolites <- distinct(metabolites) 

##Map pathways to metabolites
metabolites_pathway_df <- map_pathways_to_metabolites(metabolites)

#Save metabolites-pathways data frame
write.csv(metabolites_pathway_df,"./results/metabolites_pathway_df.csv")




##Metabolite enrichment analysis with package decoupleR

##Prepare input data
network <- metabolites_pathway_df   #input 1: pathways and their associated metabolites
network$mor <- 1                    #add new column for mode of regulation
network$likelihood <- 1             #add new column for edge likelihood

matrix <- toy_metabolomic_data         #input 2: expression matrix

##Map metabolite names to KEGG IDs
matrix$metab <- row.names(matrix)
matrix$metab <- tolower(matrix$metab)
matrix$metab <- gsub(" acid","acid",matrix$metab)
matrix$metab <- gsub(" 6-phosphate","6-phosphate",matrix$metab)
matrix$metab <- gsub("adenosyl-","adenosyl",matrix$metab)

matrix$metab[matrix$metab == "d8-valine"] <- "valine"
matrix$metab[matrix$metab == "glycerylphosphorylcholine"] <- "phosphorylcholine"
matrix$metab[matrix$metab == "stearoyl_carnitine"] <- "stearoylcarnitine"
matrix$metab[matrix$metab == "sedoheptulose_7-phosphate"] <- "sedoheptulose7-phosphate"

matrix <- merge(matrix, mapping_table, by = "metab")
row.names(matrix) <- matrix$KEGG       #Convert column to row names
matrix[c("metab", "KEGG")] <- NULL     #Delete columns metabolites and KEGG
matrix <- unique(matrix)

mat <- matrix
mat <- na.omit(mat)   
mat <- as.matrix(mat) 

##Perform enrichment analysis

#'Calculate the activity enrichment score and p-value of all pathways by using
#'the conditions in the matrix (patient samples vs metabolite expression) by
#'calculating the mean over the expression of all metabolites.

enrichment <- run_mean(mat, network,
                              .source = .data$pathway,
                              .target = .data$metabolites,
                              .mor = .data$mor, 
                              .likelihood = .data$likelihood,  
                            times = 10000,            #number of permutations
                            seed = 42,                #a single integer
                            sparse = TRUE,
                            randomize_type = "rows")  #randomize matrix

##Keep only rows with statistic "normalized_mean"
enrichment <- enrichment[enrichment$statistic != "mean", ]

##Calculate mean p-value and score for every pathway
enrichment_mean <- enrichment %>%   
                    group_by(tf) %>%    
                    summarise_at(vars(c(score, p_value)),           
                                list(name = mean))       

##Calculate adjusted p-value
enrichment_mean$adj_p_value <- p.adjust(enrichment_mean$p_value_name,
                                        method = "fdr",
                                        n = length(enrichment_mean$p_value_name))

enrichment_mean <- as.data.frame(enrichment_mean)
colnames(enrichment_mean)[colnames(enrichment_mean) == "tf"] <- "pathway" 
colnames(enrichment_mean)[colnames(enrichment_mean) == "p_value_name"] <- "p-value" 
colnames(enrichment_mean)[colnames(enrichment_mean) == "adj_p_value"] <- "Adjusted p-value" 
colnames(enrichment_mean)[colnames(enrichment_mean) == "score_name"] <- "Enrichment score" 
enrichment_mean <- enrichment_mean[ , c(1, 3, 4, 2)]  #reorder columns

##Visualize most significant pathways in a bar plot
barplot = barplot_pathways(enrichment_mean)
barplot


##Perform gene set variation analysis (GSVA)
gene_set_variation <- run_gsva(mat, network, pathway, metabolites, verbose = FALSE)

gene_set_variation_mean <- gene_set_variation %>%   
                                  group_by(tf) %>%    
                                  summarise_at(vars(score),           
                                               list(name = mean))       

gene_set_variation_mean <- gene_set_variation_mean %>% 
                                rename(pathway = tf, score = name)



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
metabolites_universe <- unique(row.names(matrix))

##Remove compartment information and "cpd:" to get pure KEGG IDs
t_table_kegg <- t_table
t_table_kegg$KEGG <- get_pure_kegg_ids(t_table_kegg$KEGG)

##Create a vector of all significant metabolites
p_value <- 2.0 #define cutoff p-value
metabolites_significant_df <- t_table_kegg[t_table_kegg$tumorVsHealthy <= -p_value | t_table_kegg$tumorVsHealthy >= p_value, ] 
metabolites_significant <- unique(metabolites_significant_df$KEGG)


##Performing the ORA with Fisher's exact test (from piano package)
sig_pathways <- runGSAhyper(genes = metabolites_significant, 
                            universe = metabolites_universe,
                            gsc = loadGSC(metabolites_pathway_df),
                            adjMethod = "fdr")
sig_pathways_df <- as.data.frame(sig_pathways$resTab)  %>% 
  tibble::rownames_to_column(var = "pathway")

##Visualize most significant pathways in a bar plot
barplot_ORA = barplot_pathways(sig_pathways_df)
barplot_ORA
