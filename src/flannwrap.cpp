// Software License Agreement (BSD License)
//
// Copyright (c) 2014, Tilo Wiklund (tilo@wiklund.co)
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//
//     Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in
//     the documentation and/or other materials provided with the
//     distribution.
//
//     The names of its contributors may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#if defined(NDEBUG)
#undef NDEBUG
#endif

//#include <boost/config.hpp>
#include "../inst/include/flanner.h"

#include <vector>
#include <memory>

//#include <RcppCommon.h>

using namespace Rcpp;
using namespace std;

//' Create a spatial index from data.frame
//'
//' Precomputes a spatial lookup structure for later nearest neighbour queries.
//'
//' @param df data.frame containing set of points (one per row)
//' @param maxleaf maximum number of nodes in a leaf, trade off between faster queries and slower build times (small maxleaf) and slower queries and faster build times (large maxleaf)
//'
//' @details nanoflann documentation suggests a value between 10 (default) and
//' 50 for query-heavy applications. As with \code{\link{lookupKNNDataFrame}}
//' this is a lower level function that performs no sanity checking (such as df
//' containing only finite numbers and no \code{NA}).
//' @export
//[[Rcpp::export]]
XPtr<FlannerIndexPair> fromDataFrame(vec2 df, const int maxleaf=10) {
    return XPtr<FlannerIndexPair>(new FlannerIndexPair(df, maxleaf));
}

//' Nearest neighbour lookup
//'
//' Perform a nearest neighbour lookup using a precomputed index structure.
//'
//' @param pair the spatial indexing structure to do lookups with respect to
//' @param points A data.frame containing the points to lookup neighbours of
//' @param k number of neighbours to find
//'
//' @details Note that this is a lower level function that does no kind of
//' sanity checking, it may fail in strange ways if \code{points} has the an
//' incorrect number of columns, and column names are not taken into account.
//' Furthermore 0 based indices are returned.
//' @export
//[[Rcpp::export]]
DataFrame lookupKNNDataFrame( XPtr<FlannerIndexPair> pair
                            , vec2 points
                            , windows_is_dumb k ) {
    DataFrameIndex& index = pair->index;
    const size_t npoints = points[0].size();
    const size_t dim = points.size();
    // ARGHH, Rcpp doesn't like long integers >_>
    // NOTE: I think nanoflann has its type parameterised over the integer type
    // to be used for indexing.
    //IntegerVector indices(npoints*k);
    vector<size_t> indices(npoints*k);
    NumericVector sqdistances(npoints*k);
    //
    // Note: This should be a unique_ptr, but I can't make
    // compileAttributes realise I want C++11, and apparently
    // CRAN doens't allow -std=c++11, so leaving it as it is for
    // now, in case I want to upload it to CRAN.
    // unique_ptr<double> point = new double[dim];
    double* point = new double[dim];
    for(size_t ipt = 0; ipt < npoints; ++ipt) {
        for(size_t i = 0; i < dim; ++i) {
            //TODO: Figure out why this is necessary
            const vector<double>& tmp = points[i];
            point[i] = tmp[ipt];
        }
        index.knnSearch( point, k
                       , &indices[k*ipt]
                       , &sqdistances[k*ipt] );
    }
    delete point;
    //
    return DataFrame::create( Named("index")=wrap(IntegerVector(indices.begin(), indices.end()))
                            , Named("distance")=wrap(sqdistances) );
}

//n[n[nRcpp::exportn]n]n
// DataFrame lookupRadius(XPtr<DataFrameIndex> index, double r) {
// }
