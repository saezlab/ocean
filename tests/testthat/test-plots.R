ocean_path <- system.file(package="ocean")

# Input
metactivity_res <- readRDS(file.path(ocean_path,
                                     "testdata/met_activity.rds"))
mean_ES_df <- metactivity_res$ES
mean_NES_df <- metactivity_res$NES

regulons_df <- readRDS(file.path(ocean_path, "testdata/met_regulons.rds"))
translated_results <- readRDS(file.path(ocean_path,
                                        "testdata/translated_res.rds"))
t_table <- readRDS(file.path(ocean_path, "testdata/meta_map.rds"))

# Test Translate Results
test_that("test_translate", {
  res1 <- translate_results(regulons_df = regulons_df,
                            t_table = t_table,
                            mapping_table = mapping_table)
  
  expect_equal(res1, translated_results)
})

# Test Single Enzyme Plots
test_that("test enzyme plots", {
  res2 <- plotMetaboliteContribution(enzyme = '4967_1738_8050_1743', 
                                      stat_df = translated_results$t_table, 
                                      metabolite_sets = translated_results$regulons_df, 
                                      contrast_index = 1, 
                                      stat_name = 't', 
                                      scaling_factor = 5, nLabels = 15)
  
  exp2 <- readRDS(file.path(ocean_path, "testdata/enzyme_plots.rds"))
  
  expect_equal(res2$scatter$data %>% arrange(contribution),
               exp2$scatter$data %>% arrange(contribution))
  
  expect_equal(res2$cumsumPlot$data %>% arrange(contribution),
               exp2$cumsumPlot$data %>% arrange(contribution))
})


test_that("test pathway plot", {
  res3 <- pathway_HM(mean_NES_df = mean_NES_df,
                     pathway_name = 'KEGG_CITRATE_CYCLE_TCA_CYCLE',
                     pathways = kegg_pathways)
  
  exp3 <- readRDS(file.path(ocean_path, "testdata/pathway_plot.rds"))
  
  expect_equal(res3$data, exp3$data)
  
})
