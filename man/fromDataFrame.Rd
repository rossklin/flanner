\name{fromDataFrame}
\alias{fromDataFrame}
\title{Create a spatial index from data.frame}
\usage{
fromDataFrame(df, maxleaf = 10L)
}
\arguments{
  \item{df}{data.frame containing set of points (one per
  row)}

  \item{maxleaf}{maximum number of nodes in a leaf, trade
  off between faster queries and slower build times (small
  maxleaf) and slower queries and faster build times (large
  maxleaf)}
}
\description{
Precomputes a spatial lookup structure for later nearest
neighbour queries.
}
\details{
nanoflann documentation suggests a value between 10
(default) and 50 for query-heavy applications. As with
\code{\link{lookupKNNDataFrame}} this is a lower level
function that performs no sanity checking (such as df
containing only finite numbers and no \code{NA}).
}

