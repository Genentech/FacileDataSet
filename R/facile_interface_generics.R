#####################################################################################################
#### The FacileAPI: Classes from the FacileVerse must implement methods on each of these generics ###
#####################################################################################################
#
##' The Facile API
##'
##' This is a stub for a manpage all about the FacileAPI
#
### General getters
#
##' @family FacileInterface
##' @export
#fetch_organism <- function(x, ...) {
#  UseMethod("fetch_organism")
#}
#
##' @export
#fetch_organism.default <- function(x, ...) {
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
##' @export
#samples.default <- function(x, ...) {
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
#
##' @family FacileInterface
##' @export
#default_assay <- function(x, ...) {
#  UseMethod("default_assay")
#}
#
##' @export
#default_assay.default <- function(x, ...) {
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
#
##' Retrieves grouping table for samples within a FacileDataSet.
##'
##' It is natural to define subgroups of samples within larger datasets.
##' This function returns grouping definitions (which we call "facets") for
##' a `FacileDataStore`.
##'
##' @family FacileInterface
##'
##' @param x An object of a class implementing the FacileInterface
##' @param name The specific facet (grouping) definition to return. Note that
##'   this parameter isn't yet used. Only one facet table was originally
##'   defined for each FacileDataSet, but we want to enable different facet
##'   definitions to be used in the future.
##' @return A `tibble` that defines the `dataset,sample_id` tuples that belong
##'   to each "facet" (group).
##' @export
#facet_frame <- function(x, name = "default", ...) {
#  UseMethod("facet_frame")
#}
#
##' @family FacileInterface
##' @export
#facet_frame.default <- function(x, name = "default", ...) {
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
#
### Filter
#
##' @family FacileInterface
##' @export
#filter_features <- function(x, ...) {
#  UseMethod("filter_features")
#}
#
##' @export
#filter_features.default <- function(x, ...) {
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
#
##' @family FacileInterface
##' @export
#filter_samples <- function(x, ..., with_covariates = FALSE) {
#  UseMethod("filter_samples")
#}
#
##' @export
#filter_samples.default <- function(x, ..., with_covariates = FALSE) {
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
#
### Fetch (<Mean Girls reference here>)
#
##' @family FacileInterface
##' @export
#fetch_samples <- function(x, samples=NULL, assay="rnaseq", ...) {
#  UseMethod("fetch_samples")
#}
#
##' @export
#fetch_samples.default <- function(x, samples = NULL, assay = "rnaseq", ...) {
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
#
##' @family FacileInterface
##' @export
#fetch_sample_statistics <- function(x, samples=NULL, semi=TRUE, assay_name='rnaseq') {
#  UseMethod("fetch_sample_statistics")
#}
#
##' @export
#fetch_sample_statistics.default <- function(x, samples=NULL, semi=TRUE, assay_name='rnaseq') {
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
#
##' @family FacileInterface
##' @export
#assay_names <- function(x, default_first=TRUE) {
#  UseMethod("assay_names")
#}
#
##' @export
#assay_names.default <- function(x, default_first = TRUE) {
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
#
##' @family FacileInterface
##' @export
#fetch_assay_data <- function(x, features, samples=NULL,
#                             assay_name=default_assay(x),
#                             normalized=FALSE, as.matrix=FALSE, ...,
#                             subset.threshold=700, aggregate.by=NULL,
#                             verbose=FALSE) {
#  UseMethod("fetch_assay_data")
#}
#
##' @export
#fetch_assay_data.default <- function(x, features, samples=NULL,
#                             assay_name=default_assay(x),
#                             normalized=FALSE, as.matrix=FALSE, ...,
#                             subset.threshold=700, aggregate.by=NULL,
#                             verbose=FALSE) {
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
#
##' @family FacileInterface
##' @export
#fetch_assay_score <- function(x, features, samples=NULL, assay_name=NULL,
#                              as.matrix=FALSE, ..., subset.threshold=700) {
#  UseMethod("fetch_assay_score")
#}
#
##' @export
#fetch_assay_score.default <- function(x, features, samples=NULL, assay_name=NULL,
#                              as.matrix=FALSE, ..., subset.threshold=700) {
#
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
#
##' @family FacileInterface
##' @export
#fetch_sample_covariates <- function(x, samples=NULL, covariates=NULL,
#                                    custom_key=Sys.getenv("USER")) {
#  UseMethod("fetch_sample_covariates")
#}
#
##' @export
#fetch_sample_covariates.default <- function(x, samples=NULL, covariates=NULL,
#                                    custom_key=Sys.getenv("USER")) {
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
#
##' @family FacileInterface
##' @export
#fetch_custom_sample_covariates <- function(x, samples=NULL, covariates=NULL,
#                                           custom_key=Sys.getenv("USER"),
#                                           file.prefix="facile") {
#  UseMethod("fetch_custom_sample_covariates")
#}
#
##' @export
#fetch_custom_sample_covariates.default <- function(x, samples=NULL, covariates=NULL,
#                                           custom_key=Sys.getenv("USER"),
#                                           file.prefix="facile") {
#  stop("The FacileAPI requires that a specific method be written for this type.")
#}
