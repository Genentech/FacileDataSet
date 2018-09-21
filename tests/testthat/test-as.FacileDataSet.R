context("as.FacileDataSet")

#es <- multiGSEA::exampleExpressionSet(do.voom=FALSE)
#esl <- list(first=es, second=es)
#colnames(esl[['second']]) <- paste0('two_', colnames(esl[['second']]))
#
#test_that("single ExpressionSet converts to FacileDataSet", {
#})
#
#test_that("list of ExpressionSets convert to FacileDataSet", {
#})

test_that("We can get pdata metadata", {
  stopifnot(requireNamespace("Biobase", quitely = TRUE))
  stopifnot(requireNamespace("survival", quietly = TRUE))
  sinfo = data.frame(a = 1:4,
                     b = survival::Surv(1:4, c(1,1,0,1)),
                     stringsAsFactors = FALSE
  )
  rownames(sinfo) = letters[1:4]
  attr(sinfo, "label") = c(a = "a is a", b = "b is b")
  vals = matrix(1:16, ncol = 4, dimnames = list(LETTERS[1:4], letters[1:4]))
  es = Biobase::ExpressionSet(vals, Biobase::AnnotatedDataFrame(sinfo))

  expect_identical(
    FacileDataSet:::pdata_metadata(es),
    list(a = list(description = "a is a"),
         b = list(description = "b is b"))
  )
})
testthat::test_that(desc = "as.FacileDataSet.R::as.FacileDataSet.ExpressionSet", code = {
  # get ExpressionSetSample
  ExpressionSetSample <-
    readRDS(
      file = system.file(
        package = "FacileDataSet",
        "extdata",
        "exampleESandSE",
        "ExpressionSetSample.RDS"
      )
    )
  # set dir for FacileDataSet
  fdsDir <- file.path(tempdir(), "fdsDir")
  unlink(fdsDir, recursive = TRUE)
  as.FacileDataSet(
    x = list(ExpressionSetSample = ExpressionSetSample),
    path = fdsDir,
    source_assay = "exprs",
    dataset_name = "ExpressionSetSample",
    assay_name = "assayDataSample",
    organism = "unspecified"
  )
  testthat::expect_identical(object = file.size(list.files(fdsDir, full.names = TRUE)),
                             expected = c(4096, 36433, 434176, 876))
  unlink(fdsDir, recursive = TRUE)
})
testthat::test_that(desc = "as.FacileDataSet.R::as.FacileDataSet.SummarizedExperiment", code = {
  # get SummarizedExperimentSample
  SummarizedExperimentSample <-
    readRDS(
      file = system.file(
        package = "FacileDataSet",
        "extdata",
        "exampleESandSE",
        "SummarizedExperimentSample.RDS"
      )
    )
  # set dir for FacileDataSet
  fdsDir <- file.path(tempdir(), "fdsDir")
  unlink(fdsDir, recursive = TRUE)
  myFds <-
    as.FacileDataSet(
      x = list(SAMPLE = SummarizedExperimentSample),
      path = fdsDir,
      dataset_name = "DEFAULT_NAME",
      assay_name = "myAssay",
      source_assay = "counts",
      assay_type = "rnaseq",
      organism = "unspecified"
    )
  testthat::expect_identical(object = file.size(list.files(fdsDir, full.names = TRUE)),
                             expected = c(4096, 36433, 405504, 520))
  unlink(fdsDir, recursive = TRUE)
})
