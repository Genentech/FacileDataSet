context("Testing conversion to Bioc Expression Containers")

FDS <- exampleFacileDataSet()
samples <- sample_covariate_tbl(FDS) %>%
  filter(variable == 'stage' & value == 'III') %>%
  select(dataset, sample_id)
genes <- local({
  out <- c("800", "1009", "1289", "50509", "2191", "2335", "5159")
  feature_info_tbl(FDS) %>%
    filter(feature_id %in% out) %>%
    collect %$%
    feature_id
})

test_that("fetch_assay_data results converted to DGEList", {
  e <- fetch_assay_data(FDS, genes, samples)
  y <- as.DGEList(e)
  expect_is(y, 'DGEList')

  ## check samples
  expect_is(y$samples, 'data.frame')
  expect_true(setequal(y$samples$sample_id, collect(samples)$sample_id))
  expect_type(y$samples$norm.factors, 'double')
  expect_type(y$samples$lib.size, 'double')
  expect_type(y$samples$dataset, 'character')
  expect_type(y$samples$sample_id, 'character')

  expect_is(y$genes, 'data.frame')
  expect_type(y$genes$feature_id, 'character')
  expect_true(setequal(y$genes$feature_id, genes))
  expect_type(y$genes$symbol, 'character')

  ## Check that counts match up in DGEList as they would from raw matrix fetch
  m <- fetch_assay_data(FDS, genes, samples, normalize=FALSE, as.matrix=TRUE)
  expect_true(setequal(rownames(m), rownames(y)))
  expect_true(setequal(colnames(m), colnames(y)))
  expect_equal(m[rownames(y), colnames(y)], y$counts)
})
