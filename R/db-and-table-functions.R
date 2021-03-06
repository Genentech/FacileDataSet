#' Query a table to identify its primary key(s)
#'
#' @export
#'
#' @param x a \code{FacileDataSet} or \code{SQLiteConnection}
#' @param table_name the name of the table to query
#' @return a character vector of primary keys
primary_key <- function(x, table_name) {
  if (is.FacileDataSet(x)) x <- x$con
  stopifnot(is(x, 'SQLiteConnection'))
  assert_string(table_name)
  info <- dbGetQuery(x, sprintf("PRAGMA table_info(%s);", table_name))
  filter(info, pk != 0)$name
}

#' Adds rows to a table in a FacileDataSet
#'
#' This function first checks the data in the target table \code{table_name}
#' to ensure that rows in \code{dat} that exist in \code{table_name} (by
#' checking the primary key) are not added.
#'
#' @export
#'
#' @param dat the \code{data.frame} of rows to add to the table, which must
#'   have a superset of columns present in the \code{table_name} that is being
#'   appended to
#' @param x the \code{FacileDataSet}
#' @param table_name the name of the table in \code{x} to add the rows of
#'   \code{dat} to.
#' @return invisibly returns the conformed version of \code{dat}.
append_facile_table <- function(dat, x, table_name) {
  stopifnot(is.FacileDataSet(x))
  target <- try(tbl(x$con, table_name), silent=TRUE)
  if (is(target, 'try-error')) stop("Unknown table to append to: ", table_name)
  dat <- conform_data_frame(dat, target)

  ## Ensure that we don't try to add existing rows into the database
  pk <- primary_key(x, table_name)
  if (length(pk)) {
    skip <- target %>%
      semi_join(dat, by=pk, copy=TRUE, auto_index=TRUE) %>%
      collect(n=Inf) %>%
      mutate(added=FALSE)
    if (nrow(skip)) {
      warning(nrow(skip), "/", nrow(dat), " features already in database",
              immediate.=TRUE)
    }
    add.me <- anti_join(dat, skip, by=pk)
    if (nrow(add.me)) {
      dbWriteTable(x$con, table_name, add.me, append=TRUE)
      add.me$added <- TRUE
    }
    dat <- bind_rows(add.me, skip)
  } else{
    dat$added <- TRUE
    dbWriteTable(x$con, table_name, dat, append=TRUE)
  }

  invisible(dat)
}

## Database Table Accessors ====================================================

#' @export
assay_info_tbl <- function(x) {
  stopifnot(is.FacileDataSet(x))
  tbl(x$con, 'assay_info') %>% set_fds(x)
}

#' @export
assay_feature_info_tbl <- function(x) {
  stopifnot(is.FacileDataSet(x))
  tbl(x$con, 'assay_feature_info') %>% set_fds(x)
}

#' @export
assay_sample_info_tbl <- function(x) {
  stopifnot(is.FacileDataSet(x))
  tbl(x$con, 'assay_sample_info') %>% set_fds(x)
}

#' @export
feature_info_tbl <- function(x, assay_name=NULL) {
  stopifnot(is.FacileDataSet(x))
  out <- tbl(x$con, 'feature_info')
  if (!is.null(assay_name)) {
    assert_string(assay_name)
    assay.info <- assay_info_tbl(x) %>%
      filter(assay == assay_name) %>%
      collect()
    if (nrow(assay.info) == 0) {
      stop("Unknown assay: ", assay_name)
    }
    afi <- assay_feature_info_tbl(x) %>%
      filter(assay == assay_name)
    out <- semi_join(out, afi, by=c('feature_type', 'feature_id'))
  }
  out %>% set_fds(x)
}

#' Mimics the old `gene_info` table.
#'
#' @export
gene_info_tbl <- function(x) {
  # TODO: This function needs to be removed and the code that relies on gene_info_tbl
  # should be updated.
  stopifnot(is.FacileDataSet(x))
  ## Columns:
  ## feature_id|feature_type|symbol|n_exons|length|source|hdf5_index
  hdf5.info <- assay_feature_info_tbl(x) %>%
    filter(assay == 'rnaseq')

  gi <- feature_info_tbl(x) %>%
    filter(feature_type == 'entrez') %>%
    select(feature_id, feature_type, symbol=name, n_exons=-1,
           length=effective_length, source) %>%
    inner_join(hdf5.info, by='feature_id') %>%
    set_fds(x)
}

#' Mimics old sample_stats table
#'
#' This function needs to be removed and the code that relies on
#' sample_stats_tbl be updated.
#' @export
sample_stats_tbl <- function(x) {
  assay_sample_info_tbl(x) %>%
    select(dataset, sample_id, libsize, normfactor) %>%
    set_fds(x)
}

#' @export
sample_covariate_tbl <- function(x) {
  stopifnot(is.FacileDataSet(x))
  tbl(x$con, 'sample_covariate') %>% set_fds(x)
}

#' @export
sample_info_tbl <- function(x) {
  stopifnot(is.FacileDataSet(x))
  tbl(x$con, 'sample_info') %>% set_fds(x)
}

#' Get/set db
#'
#' @rdname getsetdb
#' @export
#' @param x the object
#' @param db The \code{FacileDb} object
fds <- function(x) {
  if (is.FacileDataSet(x)) return(x)
  out <- attr(x, 'fds')
  if (is.null(out)) {
    warning("No FacileDataSet found in x (", class(x)[1L], ")", immediate.=TRUE)
  }
  out
}

#' @rdname getsetdb
#' @export
"fds<-" <- function(x, value) {
  UseMethod("fds<-", x)
}

#' @rdname getsetdb
#' @export
"fds<-.tbl" <- function(x, value) {
  attr(x, 'fds') <- value
  x
}

#' @rdname getsetdb
#' @export
"fds<-.data.frame" <- function(x, value) {
  attr(x, 'fds') <- value
  x
}

"fds<-.default" <- function(x, value) {
  attr(x, 'fds') <- value
  x
}

#' @rdname getsetdb
#' @export
set_fds <- function(x, value) {
  attr(x, 'fds') <- value
  x
}

## Unexported utility functions ================================================

#' Validates the bits required in a legit FacileDataSet directory.
#' @noRd
validate.facile.dirs <- function(path, data.fn, sqlite.fn, hdf5.fn, meta.fn,
                                 anno.dir) {
  if (!dir.exists(path)) {
    stop("Top level FacileData directory does not exist: ", path)
  }
  path <- normalizePath(path)
  if (!file.exists(data.fn)) {
    stop("Data file (database) does not exists", data.fn)
  } else {
    data.fn <- normalizePath(data.fn)
    if (dirname(data.fn) != path) {
      warning("Data file is not under parent directory", immediate.=TRUE)
    }
  }
  if (!file.exists(sqlite.fn)) {
    stop("Database file does not exists", sqlite.fn)
  } else {
    sqlite.fn <- normalizePath(sqlite.fn)
    if (dirname(sqlite.fn) != path) {
      warning("Database file is not under parent directory", immediate.=TRUE)
    }
  }
  if (!file.exists(hdf5.fn)) {
    warning("HDF5 file does not exists", hdf5.fn, immediate.=TRUE)
  } else {
    hdf5.fn <- normalizePath(hdf5.fn)
    if (dirname(hdf5.fn) != path) {
      warning("HDF5 file is not under parent directory", immediate.=TRUE)
    }
  }
  meta.fn <- assert_valid_meta_file(meta.fn) %>% normalizePath
  if (!dir.exists(anno.dir)) {
    stop("Directory for custom annotations does not exist: ", anno.dir)
  } else {
    anno.dir <- normalizePath(anno.dir)
    if (dirname(anno.dir) != path) {
      warning("Custom annotation directory not under parent directory.",
              immediate.=TRUE)
    }
  }

  list(path=path, data.fn=data.fn, sqlite.fn=sqlite.fn, hdf5.fn=hdf5.fn,
       meta.fn=meta.fn, anno.dir=anno.dir)
}
