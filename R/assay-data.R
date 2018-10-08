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

#' @importFrom checkmate assert_string assert_character assert_flag assert_number
#' @importFrom rhdf5 h5read
fetch_assay_data_tbl.FacileDataSet <- function(x,
                              assay_name,
                              feature_ids,
                              samples,
                              normalized = FALSE,
                              as.matrix = FALSE,
                              subset.threshold = 700,
                              aggregate.by = c("none", "ewm", "zscore"),
                              verbose = FALSE,
                              ...) {
  
  aggregate.by = match.arg(aggregate.by)
  assert_string(assay_name)
  assert_character(feature_ids, min.len=1L)
  samples <- assert_sample_subset(samples)
  assert_flag(normalized)
  assert_flag(as.matrix)
  assert_number(subset.threshold)
  
  finfo <- assay_feature_info(x, assay_name, feature_ids=feature_ids) %>%
    collect(n=Inf) %>%
    arrange(hdf5_index)
  atype <- finfo$assay_type[1L]
  ftype <- finfo$feature_type[1L]
  sinfo <- assay_sample_info(x, assay_name, samples) %>%
    mutate(samid=paste(dataset, sample_id, sep="__"))
  bad.samples <- is.na(pull(sinfo, hdf5_index))
  if (any(bad.samples)) {
    if (verbose) {
      warning(sum(bad.samples), " samples not found in `",
              assay_name, "`assay.", immediate.=TRUE)
    }
    sinfo <- sinfo[!bad.samples,,drop=FALSE]
    #sinfo <- filter(sinfo, !is.na(hdf5_index))
  }
  
  ## DEBUG: Tune chunk size?
  ## As the number of genes you are fetching increases, only subsetting
  ## out a few of them intead of first loading the whole matrix gives you little
  ## value, for instance, over all BRCA tumors (994) these are some timings for
  ## different numbers of genes:
  ##
  ##     Ngenes                   Time (seconds)
  ##     10                       0.5s
  ##     100                      0.8s
  ##     250                      1.2s
  ##     500                      2.6s
  ##     750                      6s seconds
  ##     3000                     112 seconds!
  ##     unpspecified (all 26.5k) 7 seconds!
  ##
  ## I'm using `ridx` as a hack downstream to run around the issue of slowing
  ## down the data by trying to subset many rows via hdf5 instead of loading
  ## then subsetting after (this is so weird)
  ##
  ## TODO: setup unit tests to ensure that ridx subsetting and remapping back
  ## to original genes works
  fetch.some <- nrow(finfo) < subset.threshold
  ridx <- if (fetch.some) finfo$hdf5_index else NULL
  
  dat <- sinfo %>%
    group_by(dataset) %>%
    do(res={
      ds <- .$dataset[1L]
      hd5.name <- paste('assay', assay_name, ds, sep='/')
      vals <- h5read(hdf5fn(x), hd5.name, list(ridx, .$hdf5_index))
      if (is.null(ridx)) {
        vals <- vals[finfo$hdf5_index,,drop=FALSE]
      }
      dimnames(vals) <- list(finfo$feature_id, .$samid)
      if (normalized) {
        vals <- normalize.assay.matrix(vals, finfo, ., verbose=verbose, ...)
      }
      vals
    }) %>%
    ungroup
  
  ## NOTE: We can avoid the monster matrix creation if we only want !as.matrix
  ## returns, but this makes the code easier to reason. We can come back to this
  ## to optimize for speed later. The problem is introduced when the
  ## aggregate.by parameter was introduced
  vals <- do.call(cbind, dat$res)
  
  if (nrow(vals) == 1L) {
    if (!is.null(aggregate.by) && verbose) {
      warning("No assay feature aggregation performed over single feature",
              immediate.=TRUE)
    }
    aggregate.by <- NULL
  }
  
  if (!identical(aggregate.by, "none")) {
    scores <- switch(aggregate.by,
                     ewm=eigenWeightedMean(vals, ...)$score,
                     zscore=zScore(vals, ...)$score)
    vals <- matrix(scores, nrow=1, dimnames=list('score', names(scores)))
  }
  
  if (!as.matrix) {
    vals <-
      .melt.assay.matrix(
        v = vals,
        assay_name = assay_name,
        atype = atype,
        ftype = ftype,
        finfo = finfo
      )
    if (!identical(aggregate.by,"none")) {
      vals[, feature_type := 'aggregated']
      vals[, feature_id := 'aggregated']
      vals[, feature_name := 'aggregated']
    }
    vals <- as.tbl(setDF(vals))
  }
  
  class(vals) <- c('FacileExpression', class(vals))
  set_fds(vals, x)
}

#' @import data.table
.melt.assay.matrix <- function(v, assay_name, atype, ftype, finfo) {
  v <- as.data.table(v, keep.rownames=TRUE)
  v <- melt.data.table(v, id.vars='rn', variable.factor=FALSE,
                       variable.name='sample_id')
  setnames(v, 1L, 'feature_id')
  
  v[, dataset := sub('__.*$', '', sample_id)]
  v[, sample_id := sub('^.*__', '', sample_id)]
  v[, assay := assay_name]
  v[, assay_type := atype]
  v[, feature_type := ftype]
  xref <- match(v$feature_id, finfo$feature_id)
  v[, feature_name := finfo$name[xref]]
  corder <- c('dataset', 'sample_id', 'assay', 'assay_type',
              'feature_type', 'feature_id', 'feature_name', 'value')
  setcolorder(v, corder)
  return(v)
}