#' Fetch data from single assay of choice
#'
#' @export
#' @importFrom rhdf5 h5read
#' @importFrom multiGSEA eigenWeightedMean
#' @param x A \code{FacileDataSet} object.
#' @param features a feature descriptor (data.frame with assay and feature_id
#'   columms)
#' @param samples a sample descriptor to specify which samples to return data
#'   from.
#' @param normalized return normalize or raw data values, defaults to
#'   \code{raw}
#' @param as.matrix by default, the data is returned in a long-form tbl-like
#'   result. If set to \code{TRUE}, the data is returned as a matrix.
#' @param ... parameters to pass to normalization methods
#' @param subset.threshold sometimes fetching all the genes is faster than
#'   trying to subset. We have to figure out why that is, but I've previously
#'   tested random features of different lengths, and around 700 features was
#'   the elbow.
#' @param aggregate.by do you want individual level results or geneset
#'   scores? Use 'ewm' for eigenWeightedMean, and that's all.
#' @return A lazy \code{\link[dplyr]{tbl}} object with the expression
#'   data to be \code{\link[dplyr]{collect}}ed when \code{db} is provided,
#'   otherwise a \code{tbl_df} of the results.
#' @family API
fetch_assay_data.FacileDataSet <- function(x, features = NULL, samples=NULL,
                             assay_name=default_assay(x),
                             normalized=FALSE, ..., as.matrix=FALSE,
                             subset.threshold=700, aggregate.by=c("none", "ewm", "zscore"),
                             verbose=FALSE) {
  assert_flag(as.matrix)
  assert_flag(normalized)
  assert_number(subset.threshold)
  aggregate.by = match.arg(aggregate.by)

  if (!is.null(assay_name) || is.character(features)) {
    assert_string(assay_name)
    assert_choice(assay_name, assay_names(x))
  }

  if (is.null(features)) {
    assert_string(assay_name)
    features <- assay_feature_info(x, assay_name) %>% collect(n=Inf)
  } else {
    if (is.character(features)) {
      features <- tibble(feature_id=features, assay=assay_name)
    }
    stopifnot(is(features, 'tbl') || is(features, 'data.frame'))
    if (!'assay' %in% colnames(features) || !is.character(features$assay)) {
      features$assay <- assay_name
    }
    assert_assay_feature_descriptor(features)
  }

  ## Adding check for a 0 row data.frame, because there is some chain of
  ## reactivity that fires in FacileExplorer upon FDS switching that triggers
  ## this when the previous dataset has sample filters entered. This acts
  ## as a defense to that, and also works to handle this strange case from the
  ## backend side, too -- perhaps a user will stumble on this in their analyses?
  if (is.null(samples) || (is.data.frame(samples) && nrow(samples) == 0L)) {
    samples <- fetch_samples(x)
  }
  assert_sample_subset(samples)

  assays <- unique(features$assay)
  n.assays <- length(assays)
  if (n.assays > 1L && as.matrix) {
    stop("Fetching from multiple assays requires return in melted form")
  }

  if (!identical(aggregate.by, "none")) {
    stopifnot(n.assays == 1L)
    if (!normalized) {
      warning("You probably don't want to aggregate.by on unnormalized data",
              immediate.=TRUE)
    }
  }

  out <- lapply(assays, function(a) {
    f <- filter(features, assay == a)
    .fetch_assay_data(x, a, f$feature_id, samples, normalized, as.matrix,
                      subset.threshold, aggregate.by, verbose=verbose)
  })

  if (length(out) == 1L) {
    out <- out[[1L]]
  } else if (!as.matrix) {
    ## We stop if we are asking for a matrix across multiple assays, but maybe
    ## we don't have to ... (in the future, I mean)
    out <- bind_rows(out)
  }

  out
}

.fetch_assay_data <- function(x, assay_name, feature_ids, samples,
                              normalized=FALSE, as.matrix=FALSE,
                              subset.threshold=700, aggregate.by=c("none", "ewm", "zscore"),
                              ..., verbose=FALSE) {
  #  stopifnot(is.FacileDataSet(x))
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
  bad.samples <- is.na(sinfo$hdf5_index)
  if (any(bad.samples)) {
    if (verbose) {
      warning(sum(bad.samples), " samples not found in `",
              assay_name, "`assay.", immediate.=TRUE)
    }
    sinfo <- sinfo[!bad.samples,,drop=FALSE]
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
    vals <- .melt.assay.matrix(vals, assay_name, atype, ftype, finfo)
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

.melt.assay.matrix <- function(vals, assay_name, atype, ftype, finfo) {
  vals <- as.data.table(vals, keep.rownames=TRUE)
  vals <- melt.data.table(vals, id.vars='rn', variable.factor=FALSE,
                          variable.name='sample_id')
  setnames(vals, 1L, 'feature_id')

  vals[, dataset := sub('__.*$', '', sample_id)]
  vals[, sample_id := sub('^.*__', '', sample_id)]
  vals[, assay := assay_name]
  vals[, assay_type := atype]
  vals[, feature_type := ftype]
  xref <- match(vals$feature_id, finfo$feature_id)
  vals[, feature_name := finfo$name[xref]]
  corder <- c('dataset', 'sample_id', 'assay', 'assay_type',
              'feature_type', 'feature_id', 'feature_name', 'value')
  setcolorder(vals, corder)
  vals
}

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
