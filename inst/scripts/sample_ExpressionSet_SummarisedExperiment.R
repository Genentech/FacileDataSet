# 2018.09.17
#### 1. create an instance of ExpressionSet

library(Biobase)
### assayDataSample
assayDataSample <- matrix(runif(1000), nrow=100, ncol=10)
# requred colnames
colnames(assayDataSample) <- letters[1:10]
# phenoDataSample
phenoDataSample <-
  data.frame(sex = sample(c("M", "F"), size = 10, replace = TRUE),
             type = sample(c("Control", "Case"), size = 10, replace = TRUE),
             score = round(runif(10, min = 0.6),2), stringsAsFactors = FALSE)
attr(x = phenoDataSample, which = "label") <- as.character(c("sex","type","score"))
rownames(phenoDataSample) <- letters[1:10]
### featureDataSample
featureDataSample <-
  data.frame(
    feature_type = rep("entrez", 100),
    feature_id = as.character(1:100),
    name = "character",
    meta = "character",
    effective_length = 1L,
    source = "character",
    stringsAsFactors = FALSE
  )
### ExpressionSet
esSample <-
  ExpressionSet(
    assayDataSample,
    phenoData = AnnotatedDataFrame(phenoDataSample),
    featureData = AnnotatedDataFrame(featureDataSample),
    protocolData = annotatedDataFrameFrom(assayDataSample, byrow = FALSE),
    annotation = "annotation"
  )
# !!! source_assay = "exprs"
# adata(esSample, assay = "exprs")

#### 2. create an instance of SummarisedExperiment

library(SummarizedExperiment)
nrows <- 200
ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
# assay colnames() must be NULL or equal colData rownames()
colnames(counts) <- LETTERS[1:6]
#rownames(counts) <- paste0("r",(1:dim(counts)[1]))
rownames(counts) <- sprintf("ID%03d", 1:200)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=sprintf("ID%03d", 1:200),
                     feature_type=rep("entrez", 200),
                     name="character",
                     meta="character",
                     effective_length=1L,
                     source="character"
                     )
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
S4Vectors::metadata(colData)$Treatment <- "Treatment"

seSample <- SummarizedExperiment(assays=list(counts=counts),
                     rowRanges=rowRanges, colData=colData)

#### 3. check if as.FacileDataSet() creates FacileDataSet object from se and es correctly

fdsDir <- file.path(tempdir(), "fdsDir")
myFds<- as.FacileDataSet(x = list(SAMPLE = seSample), 
                         path = fdsDir,
                         dataset_name = "DEFAULT_NAME",
                         assay_name = "myAssay", 
                         source_assay = "counts",
                         assay_type = "rnaseq", 
                         organism = "unspecified")
list.files(fdsDir)
unlink(fdsDir, recursive = TRUE)

# Warning in append_facile_table(., x, "sample_info") :
#   6/6 features already in database
# This warning is returned by as.FacileDataSet.list() -> append_facile_table()

devtools::load_all("../facileanalysis/")
devtools::load_all("../faciledataset/")
fds = FacileDataSet::exampleFacileDataSet()

