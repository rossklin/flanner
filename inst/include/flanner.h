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

#include "nanoflann.hpp"
#include <Rcpp.h>
#include <vector>

#ifdef _WIN32
typedef unsigned int windows_is_dumb;
#else
typedef size_t windows_is_dumb;
#endif

using namespace std;
using namespace Rcpp;
using namespace nanoflann;

typedef vector<vector<double> > vec2;

//TODO: Can we make DataFrame be a reference instead?
struct DataFramePointCloud {
    //TODO: Check if this gets faster by using vector<NumericVector>
    // I suspect it would, unless DataFrame does some magic to avoid
    // runtime type checking.
    vec2 points;
    size_t ndim;
    size_t count;

    DataFramePointCloud(vec2 points_) :
        points(points_), ndim(points_.size()), count(points_[0].size()) {
    }

    inline size_t kdtree_get_point_count() const { return count; }

    inline double kdtree_distance( const double *query
                                 , const size_t index
                                 , size_t size ) const {
        double accum = 0;
        for(size_t dim = 0; dim < ndim; ++dim) {
            //TODO: Try to figure out why this is required
            const vector<double>& tmp = points[dim];
            double d = tmp[index] - query[dim];
            accum += d*d;
        }
        return accum;
    }

    inline double kdtree_get_pt(const size_t index, int dim) const {
        //TODO: Try to figure out why this is required
        const vector<double>& tmp = points[dim];
        return tmp[index];
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const { return false; }
};

typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<double, DataFramePointCloud>,
    DataFramePointCloud> DataFrameIndex;

// Note, required to get sane memory management with XPtr
struct FlannerIndexPair {
    DataFramePointCloud cloud;
    DataFrameIndex index;
    // TODO: Can we somehow take df by reference?
    FlannerIndexPair(vec2 df, const int maxleaf)
        : cloud(df)
        , index(df.size(), cloud, KDTreeSingleIndexAdaptorParams(maxleaf)){
        index.buildIndex();
    }
};
