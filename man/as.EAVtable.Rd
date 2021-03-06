% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entity-attribute-value.R
\name{as.EAVtable}
\alias{as.EAVtable}
\title{Convert a \code{pData} data.frame to a melted EAV table}
\usage{
as.EAVtable(x, eav_metadata = NULL, covariate_def = list())
}
\arguments{
\item{x}{a wide \code{pData} data.frame}

\item{eav_metadata}{the list-of-list covariate encodings for the EAV table
of the type that is generated by \code{\link[=eav_metadata_create]{eav_metadata_create()}}}

\item{covariate_def}{passed to \code{\link[=eav_metadata_create]{eav_metadata_create()}} if \code{eav_metadata}
parameter is not provided.}
}
\value{
a melted EAV table from \code{x}
}
\description{
Transforms a wide \code{pData} data.frame into a melted EAV table for use in
a \code{FacileDataSet}. This function requires the list-of-list encodings that
are generated from \code{\link[=eav_metadata_create]{eav_metadata_create()}} to do its thing. The caller can
provide their own encoding via the \code{eav_metadata} parameter, otherwise a
default one will be generated.
}
