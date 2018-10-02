#' @export
assay_names.FacileDataSet <- function(x, default_first=TRUE) {
  anames <- assay_info_tbl(x) %>% collect %$% assay
  if (default_first && length(anames) > 1L) {
    dassay <- default_assay(x)
    anames <- intersect(c(dassay, setdiff(anames, dassay)), anames)
  }
  anames
}

## helper function to fetch_assay_data
normalize.assay.matrix <- function(vals, feature.info, sample.info,
                                   log=TRUE, prior.count=5, ...,
                                   verbose=FALSE) {
  stopifnot(
    nrow(vals) == nrow(feature.info),
    all(rownames(vals) == feature.info$feature_id),
    ncol(vals) == nrow(sample.info),
    all(colnames(vals) == sample.info$samid),
    is.character(feature.info$assay_type),
    length(unique(feature.info$assay_type)) == 1L,
    is.numeric(sample.info$libsize), is.numeric(sample.info$normfactor))
  atype <- feature.info$assay_type[1L]
  libsize <- sample.info$libsize * sample.info$normfactor
  if (atype == 'rnaseq') {
    # we assume these are units that are at the count level
    out <- edgeR::cpm(vals, libsize, log=log, prior.count=prior.count)
  } else if (atype == "tpm") {
    # someone processed their data with salmon or kallisto and wanted to store
    # tpm. Normalizing this is just log2(val + prior.count)
    out <- log2(vals + prior.count)
  } else {
    if (verbose) {
      warning("No normalization procedure for ", atype, " assay",
              immediate.=TRUE)
    }
    out <- vals
  }
  out
}
