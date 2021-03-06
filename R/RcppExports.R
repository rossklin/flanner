# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Create a spatial index from data.frame
#'
#' Precomputes a spatial lookup structure for later nearest neighbour queries.
#'
#' @param df data.frame containing set of points (one per row)
#' @param maxleaf maximum number of nodes in a leaf, trade off between faster queries and slower build times (small maxleaf) and slower queries and faster build times (large maxleaf)
#'
#' @details nanoflann documentation suggests a value between 10 (default) and
#' 50 for query-heavy applications. As with \code{\link{lookupKNNDataFrame}}
#' this is a lower level function that performs no sanity checking (such as df
#' containing only finite numbers and no \code{NA}).
#' @export
fromDataFrame <- function(df, maxleaf = 10L) {
    .Call('flanner_fromDataFrame', PACKAGE = 'flanner', df, maxleaf)
}

#' Nearest neighbour lookup
#'
#' Perform a nearest neighbour lookup using a precomputed index structure.
#'
#' @param pair the spatial indexing structure to do lookups with respect to
#' @param points A data.frame containing the points to lookup neighbours of
#' @param k number of neighbours to find
#'
#' @details Note that this is a lower level function that does no kind of
#' sanity checking, it may fail in strange ways if \code{points} has the an
#' incorrect number of columns, and column names are not taken into account.
#' Furthermore 0 based indices are returned.
#' @export
lookupKNNDataFrame <- function(pair, points, k) {
    .Call('flanner_lookupKNNDataFrame', PACKAGE = 'flanner', pair, points, k)
}

