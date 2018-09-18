#### 2018.09.18
#### 1. Set environment and load FacileDataSet example

library(FacileAnalysis)
library(FacileDataSet)

# load sample FacileDataSet
x = FacileDataSet::exampleFacileDataSet()
# create input data.frame for as.ExpressionSet and as.SummarizedExperiment
.fds=FacileAnalysis::fds(x)
x <- FacileAnalysis::samples(x) %>% collect(n=Inf) %>% FacileAnalysis::set_fds(.fds)

#### 2. Create ExpressionSet from FacileDataSet

library("Biobase")

### create ExpressionSet from data.frame
esSampleFE <- as.ExpressionSet.data.frame(x)

### prepare input for sample ExpressionSet
# SampleSize
SamSize <- 1000
## assayData:
assayDataSample <- adata(esSampleFE)[1:SamSize,]
## phenoData
# 'dataset' and 'sample_id' columns have to be created with as.FacileDataSet
phenoDataSample <- dplyr::select(pdata(esSampleFE), -dataset, -sample_id)
# add 'label' attribute
attr(x = phenoDataSample, which = "label") <- as.character(colnames(phenoDataSample))
## featureData
featureDataSample <- fdata(esSampleFE)[1:SamSize,]
# adding missing columns
featureDataSample$name=featureDataSample$seqnames
featureDataSample$effective_length <- featureDataSample$length
## create new ExpressionSet object
ExpressionSetSample <-
  ExpressionSet(
    assayData = assayDataSample,
    phenoData = AnnotatedDataFrame(phenoDataSample),
    featureData = AnnotatedDataFrame(featureDataSample),
    protocolData = annotatedDataFrameFrom(assayDataSample, byrow = FALSE),
    annotation = "annotation"
  )
## save sample ExpressionSetSample
saveRDS(
  object = ExpressionSetSample,
  file = file.path("inst", "extdata", "exampleESandSE", "ExpressionSetSample.RDS")
)

#### 3. Create ExpressionSet from FacileDataSet

library("SummarizedExperiment")

### create SummarizedExperiment from data.frame
seSampleFE <- FacileDataSet::as.SummarizedExperiment.data.frame(x)
### prepare input for sample SummarizedExperiment
# SampleSize
SamSize <- 1000
## assays
assaysSample <- assay(seSampleFE)[1:SamSize,]
## rowRanges
rowDataSample <- SummarizedExperiment::rowData(seSampleFE)[1:SamSize,]
rowDataSample$name <- rowDataSample$seqnames
rowDataSample$effective_length <- rowDataSample$length
## DataSample
colDataSample <- SummarizedExperiment::colData(seSampleFE)
colDataSample <- colDataSample[,!colnames(colDataSample) %in% c("dataset", "sample_id")]
S4Vectors::metadata(colDataSample) <- mapply(colnames(colDataSample), FUN =list)
### create new SummarizedExperiment object
SummarizedExperimentSample <- SummarizedExperiment(
  assays = list(counts = assaysSample),
  rowData = rowDataSample,
  colData = colDataSample
)
## save sample ExpressionSetSample
saveRDS(
  object = SummarizedExperimentSample,
  file = file.path("inst", "extdata", "exampleESandSE", "SummarizedExperimentSample.RDS")
)
