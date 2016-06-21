## Software License Agreement (BSD License)
##
## Copyright (c) 2014, Tilo Wiklund (tilo@wiklund.co)
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##     Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##     Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##
##     The names of its contributors may not be used to endorse or promote products
##     derived from this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#' Compute lookup index
#'
#' Precompute a nearest neighbours lookup index for later queries.
#'
#' @param df data.frame containing points
#' @param cols columns with respect to which to do lookups (defaults to all columns)
#' @param may.contain.NA whether the data.frame can contain incomplete rows (which will be ignored)
#' @param maxleaf maximum number of items in each leaf node (nanoflann suggests between 10 and 50, defaulting to 10)
#'
#' @details Returns \code{df} with additional class "flanner" and any attributes
#' required to make efficient nearest neighbour lookups using
#' \code{\link{knn_lookup_rows}} or \code{\link{knn_lookup}}.
#'
#' Any rows containing NA are ignored, but row numbers produced by lookup
#' functions will take NA values into account when computing indices (that is to
#' say they will return row numbers of \code{df}, *not* row numbers of \code{df}
#' after incomplete rows have been removed).
#'
#' @export
flanner <- function( df, cols=NULL
                   , may.contain.NA=TRUE
                   , maxleaf=10 ) {
    if(inherits(df, "flanner")) {
        warning("df already has an index.")
        return(df)
    }
    #
    df2 <- if(!is.null(cols)) as.data.table(df)[,cols,with=FALSE] else df
    cols <- maybe(cols, colnames(df2))
    if(length(cols) < 1) stop("Spatial index must be built with respect to at least one column!")
    ## If df contains NA values (and we allow such) create a translation
    ## table translating from row numbers of the NA-free data.frame to those
    ## ofthe original data.frame. NULL indicates no such lookup is required.
    index.lookup.table <- if(may.contain.NA) {
        is.complete <- complete.cases(df2)
        if(all(is.complete)) NULL
        else which(is.complete)
    }
    #
    spatial.index <- if(!is.null(index.lookup.table)) {
        fromDataFrame(as.data.table(df2)[complete.cases,with=FALSE], maxleaf)
    } else {
        fromDataFrame(df2, maxleaf)
    }
    #
    ## result <- structure( df
    ##                    , flanner.lookup.table=index.lookup.table
    ##                    , flanner.spatial.index=spatial.index
    ##                    , flanner.columns=cols )
    result <- copy(df)
    ### TODO: I think class can be set with setattr, but not sure, whence the extra copy...
    ### class(result) <- c("flanner", class(df))
    ### NOTE: Needed in case of data.table
    ### result <- copy(result)
    setattr(result, "class", c("flanner", class(df)))
    setattr(result, "flanner.lookup.table", index.lookup.table)
    setattr(result, "flanner.spatial.index", spatial.index)
    setattr(result, "flanner.columns", cols)
    result
}

#' Nearest neighbour rows
#'
#' Lookup row numbers of nearest neighbours in a data.frame.
#'
#' @param df data.frame, preferably indexed using \code{flanner}, in which to look for neighbours
#' @param points data.frame containing coordinates to lookup neighbours of
#' @param k number of neighbours to find
#' @param ignore.colnames whether to interpret \code{points} by column numbers of column names (defaults to the latter)
#' @param square.distance whether to return distances squared (default, which is faster)
#'
#' @details the returned row number vector has a "distance" attribute containing
#' the (possibly squared) euclidean (L^2) distances.
#'
#' @export
knn_lookup_rows <- function( df, points, k
                           , ignore.colnames=FALSE
                           , square.distance=TRUE, ... ) {
    if(!inherits(df, "flanner")) {
        warning("df is does not have a precomputed index, will compute one. Consider precomputing one using 'flanner' if you're running multiple lookups.")
        df <- flanner(df, ...)
    }
    lookup.table <- attr(df, "flanner.lookup.table")
    spatial.index <- attr(df, "flanner.spatial.index")
    #
    points2 <- if(ignore.colnames) points else as.data.table(points)[,attr(df, "flanner.columns"),with=FALSE]
    neighbours <- lookupKNNDataFrame(spatial.index, points2, k)
    #
    if(!is.null(lookup.table)) {
        neighbours$index <- lookup.table[neighbours$index+1]
    } else {
        neighbours$index <- neighbours$index+1
    }
    if(!square.distance) neighbours$distance <- sqrt(neighbours$distance)
    #
    structure(neighbours$index, distance=neighbours$distance)
}

#' Find nearest neighbour
#'
#' Lookup nearest neighbours data in a data.frame. Modified by Ross
#' Linscott so as.data.table does not destroy class attributes of df.
#'
#' @param df data.frame, preferably indexed using \code{flanner}, in which to look for neighbours
#' @param points data.frame containing coordinates to lookup neighbours of
#' @param k number of neighbours to find
#' @param ignore.colnames whether to interpret \code{points} by column numbers of column names (defaults to the latter)
#' @param square.distance whether to return distances squared (default, which is faster)
#' @param df.cols columns from df to return (defaults to all)
#' @param points.cols columns from points to return (defaults to any not used for lookup) \code{c()} for none
#' @param distance.name column name containing distances in result, \code{NULL} to not include
#'
#' @details The returned data frame contains all columns specified in df.cols and pointer.cols
#' as well as a column with distances if \code{distance.name} is non-null.
#'
#' @export
knn_lookup <- function( df, points, k
                      , ignore.colnames=FALSE
                      , square.distance=TRUE
                      , df.cols = NULL
                      , points.cols = NULL
                      , distance.name="distance" ) {
    columns <- maybe(attr(df, "flanner.columns"), colnames(df))
    #
    neighbours <- knn_lookup_rows(df, points, k, ignore.colnames, square.distance)
    #
    keep.df.cols <- maybe(df.cols, colnames(df))
    keep.points.cols <- maybe(points.cols, if(!ignore.colnames) {
        setdiff(colnames(points), union(columns, keep.df.cols))
    } else {
        colnames(points)[-seq_along(columns)]
    })
    #
    df2 <- as.data.table(copy(df))[neighbours, keep.df.cols, with=F]
    #
    cols.copy <- if(length(keep.points.cols)>0) rep(seq_len(nrow(points)), each=k)
    for(col in keep.points.cols) {
        # For some reason this is quite a bit faster than a cbind
        df2[[col]] <- points[[col]][cols.copy]
    }
    #
    if(!is.null(distance.name)) {
        df2[[distance.name]] <- attr(neighbours, "distance")
    }
    #
    df2
}
