\name{lookupKNNDataFrame}
\alias{lookupKNNDataFrame}
\title{Nearest neighbour lookup}
\usage{
lookupKNNDataFrame(pair, points, k)
}
\arguments{
  \item{pair}{the spatial indexing structure to do lookups
  with respect to}

  \item{points}{A data.frame containing the points to
  lookup neighbours of}

  \item{k}{number of neighbours to find}
}
\description{
Perform a nearest neighbour lookup using a precomputed
index structure.
}
\details{
Note that this is a lower level function that does no kind
of sanity checking, it may fail in strange ways if
\code{points} has the an incorrect number of columns, and
column names are not taken into account. Furthermore 0
based indices are returned.
}

