\name{knn_lookup}
\alias{knn_lookup}
\title{Find nearest neighbour}
\usage{
knn_lookup(df, points, k, ignore.colnames = FALSE, square.distance = TRUE,
  df.cols = NULL, points.cols = NULL, distance.name = "distance")
}
\arguments{
  \item{df}{data.frame, preferably indexed using
  \code{flanner}, in which to look for neighbours}

  \item{points}{data.frame containing coordinates to lookup
  neighbours of}

  \item{k}{number of neighbours to find}

  \item{ignore.colnames}{whether to interpret \code{points}
  by column numbers of column names (defaults to the
  latter)}

  \item{square.distance}{whether to return distances
  squared (default, which is faster)}

  \item{df.cols}{columns from df to return (defaults to
  all)}

  \item{points.cols}{columns from points to return
  (defaults to any not used for lookup) \code{c()} for
  none}

  \item{distance.name}{column name containing distances in
  result, \code{NULL} to not include}
}
\description{
Lookup nearest neighbours data in a data.frame.
}
\details{
The returned data frame contains all columns specified in
df.cols and pointer.cols as well as a column with distances
if \code{distance.name} is non-null.
}

