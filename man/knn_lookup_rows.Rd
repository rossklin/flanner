\name{knn_lookup_rows}
\alias{knn_lookup_rows}
\title{Nearest neighbour rows}
\usage{
knn_lookup_rows(df, points, k, ignore.colnames = FALSE,
  square.distance = TRUE, ...)
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
}
\description{
Lookup row numbers of nearest neighbours in a data.frame.
}
\details{
the returned row number vector has a "distance" attribute
containing the (possibly squared) euclidean (L^2)
distances.
}

