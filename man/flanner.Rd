\docType{package}
\name{flanner}
\alias{flanner}
\alias{flanner-package}
\title{flanner.}
\usage{
flanner(df, cols = NULL, may.contain.NA = TRUE, maxleaf = 10)
}
\arguments{
  \item{df}{data.frame containing points}

  \item{cols}{columns with respect to which to do lookups
  (defaults to all columns)}

  \item{may.contain.NA}{whether the data.frame can contain
  incomplete rows (which will be ignored)}

  \item{maxleaf}{maximum number of items in each leaf node
  (nanoflann suggests between 10 and 50, defaulting to 10)}
}
\description{
flanner.

Precompute a nearest neighbours lookup index for later
queries.
}
\details{
Returns \code{df} with additional class "flanner" and any
attributes required to make efficient nearest neighbour
lookups using \code{\link{knn_lookup_rows}} or
\code{\link{knn_lookup}}.

Any rows containing NA are ignored, but row numbers
produced by lookup functions will take NA values into
account when computing indices (that is to say they will
return row numbers of \code{df}, *not* row numbers of
\code{df} after incomplete rows have been removed).
}

