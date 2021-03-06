% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_FDS_pkg.R
\name{create_FDS_pkg}
\alias{create_FDS_pkg}
\title{Create a Facile___DataSet Package from a list of SummarizedExperiments}
\usage{
create_FDS_pkg(data_list = NULL, slug = NULL, version = "0.0.1",
  parent_path = ".", covariates = NULL, cov_metadata = NULL,
  data_metadata = NULL, assay_name = "rnaseq",
  source_assay = "counts", organism = "Homo sapiens")
}
\arguments{
\item{data_list}{named \code{list} of data \code{SummarizedExperiment} objects to be converted into a \code{Facile___DataSet} package.}

\item{slug}{name to use in the name \code{Facile___DataSet} package.}

\item{version}{version to label the new \code{Facile___DataSet} package.}

\item{parent_path}{directory in which the package will be created. Must already exist.}

\item{covariates}{covariates to extract from the listed data structures. Must exist in each dataset.}

\item{cov_metadata}{nested \code{list} of covariate metadata. Must contain at least label, description, type and class for each covariate.}

\item{data_metadata}{named `list`` of metadata for the datasets. Must contain at least url and description for each dataset.}

\item{assay_name}{label to be used within FacileExplorer. Default "rnaseq".}

\item{organism}{label to be used within FacileExplorer. Default "Homo Sapiens".}
}
\value{
Builds and installs a \code{Facile___DataSet} package, leaving a source tarball in \code{parent_path} directory.
}
\description{
Accepts a list of datasets (SummarizedExperiments) and metadata and
builds/installs a Facile___DataSet based on this.
}
\examples{
\dontrun{
ngs114 <- ep.project.from.id("PRJ0013166", standard.checks = FALSE)
ngs171 <- ep.project.from.id("PRJ0013155", standard.checks = FALSE)

data_list <- list(
  NGS114 = ep.SummarizedExperiment(ngs114, attach.annot = TRUE),
  NGS171 = ep.SummarizedExperiment(ngs171, attach.annot = TRUE)
)

cov_metadata <- list(
  primary_tissue = list(label = "Primary Tissue",
                        description = "Primary tissue source",
                        type = "tumor_classification",
                        class = "categorical"
  ),
  gender = list(label = "Gender/Sex",
                description = "In the 'ratio between chrX:chrY' sense.",
                type = "clinical",
                class = "categorical"
  ),
  diagnosis = list(label = "Tissue Diagnosis",
                   description = "Tissue metaclass oncology",
                   type = "tumor_classification",
                   class = "categorical"
  ),
  ethnicity = list(label = "Ethnicity",
                   description = "Patient ethnicity",
                   type = "clinical",
                   class = "categorical"
  ),
  tissue = list(label = "Tissue Group",
                description = "Rollup of tissue type to defined vocab",
                type = "tumor_classification",
                class = "categorical"
  )
)

data_metadata <- list(
  NGS114 = list(url = "http://gene.com", description = "This is NGS114"),
  NGS171 = list(url = "http://gene.com", description = "This is NGS171")
)

create_FDS_pkg(data_list = data_list,
               slug = "GCell",
               version = "0.0.1",
               parent_path = "~/FacileVerse",
               covariates = c("PRIMARY_TISSUE", "GENDER", 
                              "TISSUE_DIAGNOSIS", "ETHNICITY", 
                              "TISSUE_METACLASS_ONCOLOGY"),
               cov_metadata = cov_metadata,
               data_metadata = data_metadata)

}
}
