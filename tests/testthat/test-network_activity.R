ocean_path <- system.file(package="ocean")

# TCA_network_nocofact equivalent
filt_network <- readRDS(file.path(ocean_path,  
                                  "testdata/network_filt.rds"))
t_table <- readRDS(file.path(ocean_path, "testdata/meta_map.rds"))

penalty_min <- 6 #minimum 1 and integer
penalty_max <- 6 #maximum 9 and integer

test_that("test_network", {
  res1 <- 
    model_to_pathway_sif(pathway_to_keep = c("Citric acid cycle",
                                             "Glycolysis/gluconeogenesis",
                                             "Pyruvate metabolism",
                                             "Transport, mitochondrial")) %>%
    remove_cofactors() %>%
    compress_transporters() %>%
    split_transaminases()
  
  expect_equal(res1, filt_network)
    
})


test_that("test_forest", {
  
  # Get enzymes from the filtered network
  enzymes <- unique(filt_network$attributes$V1)
  enzymes <- enzymes[!grepl("_[clxmenr]$",enzymes)]
  
  TCA_forest <- forestMaker(enzymes,
                            filt_network$reaction_network,
                            branch_length = c(3,3))
  
  res2 <- # equivalent to regulons_df
    prepare_metabolite_set(penalty_range = penalty_min:penalty_max,  
                           forest = TCA_forest,
                           measured_metabolites = t_table$KEGG) %>%
    condense_metabolite_set() %>%
    prepare_regulon_df(., 6, c(0,1)) %>%
    arrange(weight, targets, set)
  
  exp2 <- readRDS(file.path(ocean_path, "testdata/met_regulons.rds")) %>%
    arrange(weight, targets, set)
  
  expect_equal(res2, exp2)

})


test_that("test_activity", {
  
  res3 <- metactivity(metabolomic_t_table = t_table, 
                      regulons_df = readRDS(file.path(ocean_path, "testdata/met_regulons.rds")), 
                      compartment_pattern = "_[a-z]$", 
                      k = 1000)
  
  exp3 <- readRDS(file.path(ocean_path, "testdata/met_activity.rds"))

  
  expect_equal(exp3$NES %>% arrange(row.names(.)),
               res3$NES %>% arrange(row.names(.)))
})
