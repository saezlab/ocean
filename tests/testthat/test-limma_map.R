ocean_path <- system.file(package="ocean")

unique(toy_targets$condition)
comparisons <- list('tumorVsHealthy' = c(1,-2))
limmaRes <- readRDS(file.path(ocean_path, "testdata/limma_res.rds"))
t_table <- readRDS(file.path(ocean_path, "testdata/ttop_res.rds"))

# Test Limma
test_that("test_limma", {
  res1 <- runLimma(measurements = toy_metabolomic_data,
                   targets = toy_targets,
                   comparisons = comparisons)
  
  expect_equal(res1, limmaRes)
})

# Test ttop formatter
test_that("test_ttop", {
  res2 <- ttop_list_to_t_table(
    limma_res_to_ttop_list(limma_res = limmaRes,
                           comp_names = names(comparisons),
                           number = length(toy_metabolomic_data[,1]),
                           adjust.method = "fdr"))
  expect_equal(res2, t_table)
  
})


# Test metabolic mappper
test_that("meta_map", {
  exp3 <- readRDS(file.path(ocean_path, "testdata/meta_map.rds"))
  
  res3 <- t_table_metactivity_input_formater(metabolomic_t_table = t_table,
                                             mapping_table = mapping_table,
                                             affixes = c("c","l","x","m","e","n","r"))
  expect_equal(res3, exp3)
  
})


