% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FacileDataSet.R
\name{FacileDataSet}
\alias{FacileDataSet}
\title{Instantiates a FacileDataSet object from disk.}
\usage{
FacileDataSet(path, data.fn = file.path(path, "data.sqlite"),
  sqlite.fn = file.path(path, "data.sqlite"), hdf5.fn = file.path(path,
  "data.h5"), meta.fn = file.path(path, "meta.yaml"),
  anno.dir = file.path(path, "custom-annotation"), cache_size = 80000,
  db.loc = c("reference", "temporary", "memory"), ...)
}
\arguments{
\item{path}{The path to the FacileData repository}

\item{data.fn}{A custom path to the database (probably don't mess with this)}

\item{sqlite.fn}{name of SQLite data file in FacileDataSet}

\item{hdf5.fn}{name of HDF5 data file in FacileDataSet}

\item{meta.fn}{name of metadata YAML data file in FacileDataSet}

\item{anno.dir}{A directory to house custom annotations/sample covariates}

\item{cache_size}{A custom paramter for the SQLite database}

\item{db.loc}{single character, location for the data}

\item{...}{other args to pass down, not used at the moment}

\item{covdef.fn}{A custom path to the yaml file that has covariate mapping info}
}
\value{
a `FacileDataSet` object
}
\description{
The `FacileDataSet` is a reference data storage implementation that
implements the **FacileData Access API**. It facilitates the storage and
retrieval of large amounts of data by leveraging a SQLite database to store
sample- and feature-level metadata ("`pData`" and "`fData`"), and an HDF5
file to store all of the dense assay (matrix) data (gene counts, microarray
intensities, etc.) run over the samples.
}
\details{
A `FacileDataSet` is materialized on disk by a well-structured directory,
which minimally includes the following items:

1. A `data.sqlite` SQLite database that stores feature and sample metadata
2. A `data.h5` HDF5 file that stores a multitude of dense assay matrices that
   are generated from the assays performed on the samples in the
   `FacileDataSet`.
3. A `meta.yaml` file tha contains informaiton about the `FacileDataSet`.
   To better understand the structure and contents of this file, you can
   refer to the following:
    a. The included `testdata/expected-meta.yaml` file, which is an
       exemplar file for [exampleFacileDataSet()].
    b. The help file provided by the [eav_metadata_create()] function, which
       describes in greater detail how we track a dataset's sample-level
       covariates (aka, "pData" in the bioconductor world).
   In the meantime, a short description of the entries found in the
   `meta.yaml` file is provded here:
    - `name`: the name of the dataset (ie. `"FacileTCGADataSet"`)
    - `organism`: `"Homo sapiens"`, `"Mus musculus"`, ec.
    - `default_assay`: the name of the assay to use by default if none is
      specified in calls to [fetch_assay_data()], [with_assay_data()], etc.
      (kind of like how `"exprs"` is the default assay used when working with
      a [Biobase::ExpressionSet])
    - `datasets`: a section tha enumerates the datases included internally.
      The datasets are further enumerated.
    - `sample_covariates`: a section that enumerates the covariatets that
      are tracked over the samples inside the `FacileDataSet` (ie. a mapping
      of the `pData` for the samples). Reference [eav_metadata_create()]
      for more information.
4. A `custom-annotation` directory, which stores custom `sample_covariate`
   (aka "pData") informaiton that analysts can identify and describe during
   the course of an analysis, or even add from external sources. Although
   this directory is required in the directory structure of a valid
   `FacileDataSet`, the `FacileDataSet()` constructor can be called with
   a custom `anno.dir` parameter so that custom annotations are stored
   elsewhere.
}
\examples{
fn <- system.file("extdata", "exampleFacileDataSet", package = "FacileDataSet")
fds <- FacileDataSet(fn)
}
\seealso{
Other FacileDataSet: \code{\link{dbfn}},
  \code{\link{hdf5fn}}, \code{\link{meta_file}}
}
\concept{FacileDataSet}
