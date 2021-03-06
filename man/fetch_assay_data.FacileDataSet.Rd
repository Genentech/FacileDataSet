% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assay-data.R
\name{fetch_assay_data.FacileDataSet}
\alias{fetch_assay_data.FacileDataSet}
\title{Fetch data from single assay of choice}
\usage{
\method{fetch_assay_data}{FacileDataSet}(x, features, samples = NULL,
  assay_name = default_assay(x), normalized = FALSE,
  as.matrix = FALSE, ..., subset.threshold = 700,
  aggregate.by = NULL, verbose = FALSE)
}
\arguments{
\item{x}{A \code{FacileDataSet} object.}

\item{features}{a feature descriptor (data.frame with assay and feature_id
columms)}

\item{samples}{a sample descriptor to specify which samples to return data
from.}

\item{assay_name}{the name of the assay}

\item{normalized}{return normalize or raw data values, defaults to
\code{raw}}

\item{as.matrix}{by default, the data is returned in a long-form tbl-like
result. If set to \code{TRUE}, the data is returned as a matrix.}

\item{...}{parameters to pass to normalization methods}

\item{subset.threshold}{sometimes fetching all the genes is faster than
trying to subset. We have to figure out why that is, but I've previously
tested random features of different lengths, and around 700 features was
the elbow.}

\item{aggregate.by}{do you want individual level results or geneset
scores? Use 'ewm' for eigenWeightedMean, and that's all.}
}
\value{
A lazy \code{\link[dplyr]{tbl}} object with the expression
  data to be \code{\link[dplyr]{collect}}ed when \code{db} is provided,
  otherwise a \code{tbl_df} of the results.
}
\description{
Fetch data from single assay of choice
}
\seealso{
Other API: \code{\link{fetch_assay_score.FacileDataSet}},
  \code{\link{fetch_custom_sample_covariates.FacileDataSet}},
  \code{\link{fetch_organism.FacileDataSet}},
  \code{\link{fetch_sample_covariates.FacileDataSet}},
  \code{\link{fetch_sample_statistics.FacileDataSet}},
  \code{\link{fetch_samples.FacileDataSet}},
  \code{\link{filter_features.FacileDataSet}},
  \code{\link{filter_samples.FacileDataSet}},
  \code{\link{samples.FacileDataSet}}
}
\concept{API}
