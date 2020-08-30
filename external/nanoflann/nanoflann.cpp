#include <nanoflann.hpp>
#include "utils.h"

using namespace std;
using namespace nanoflann;

template <typename num_t>
void kdtree(const size_t N, num_t *points)
{

  PointCloud<num_t> cloud;

  cloud.pts.resize(N);
  for (size_t i = 0; i < N; i++) {
    cloud.pts[i].x = points[i][0];
    cloud.pts[i].y = points[i][1];
    cloud.pts[i].z = points[i][2];
  }

  // construct a kd-tree index:
  typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<num_t, PointCloud<num_t> > ,
    PointCloud<num_t>,
    3 /* dim */
    > my_kd_tree_t;

  my_kd_tree_t index(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf, to be tuned */) );

  index.buildIndex();

  const num_t search_radius = static_cast<num_t>(0.1);
  std::vector<std::pair<size_t,num_t> >   ret_matches;

  nanoflann::SearchParams params;

  size_t maxNeigh = 0;
  for (size_t i = 0; i < N; i++) {

    const num_t query_pt[3] = { cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z};

    const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);
    maxNeigh = maxNeigh > nMatches ? maxNeigh : nMatches;

  }

}
