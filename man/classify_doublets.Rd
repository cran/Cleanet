% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{classify_doublets}
\alias{classify_doublets}
\title{Classify doublets (or multiplets) based on component singlets.}
\usage{
classify_doublets(cleanet_res, singlet_clas, max_multi = 4)
}
\arguments{
\item{cleanet_res}{The output of a call to the cleanet function.}

\item{singlet_clas}{An array giving a classification of the singlets, whose
length must match the number of singlet events returned in cleanet_res.}

\item{max_multi}{The highest cardinality of a multiplet to be considered.}
}
\value{
An array with the same length as the number of doublets found
in cleanet_res, specifying the composition of each doublet.
}
\description{
Extends a classification of singlets into a classification
of doublets.
}
\examples{
path <- system.file("extdata", "df_mdipa.csv", package="Cleanet")
df_mdipa <- read.csv(path, check.names=FALSE)
cols <- c("CD45", "CD123", "CD19", "CD11c", "CD16",
          "CD56", "CD294", "CD14", "CD3", "CD20",
          "CD66b", "CD38", "HLA-DR", "CD45RA",
          "DNA1", "DNA2")
cleanet_res <- cleanet(df_mdipa, cols, cofactor=5)
singlet_clas <- df_mdipa$label[which(cleanet_res$status!="Doublet")]
doublet_clas <- classify_doublets(cleanet_res, singlet_clas)
}
