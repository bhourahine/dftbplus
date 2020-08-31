#include <nanoflann.hpp>
#include "utils.h"
#include <iostream>

using namespace std;
using namespace nanoflann;

extern "C" {

  void kdtree(int *N, double *search_radius, double* points)
  {

    PointCloud<double> cloud;

    cloud.pts.resize(*N);

    for (int i = 0; i < *N; i++) {
      cloud.pts[i].x = points[3*i];
      cloud.pts[i].y = points[3*i+1];
      cloud.pts[i].z = points[3*i+2];
    }

    // construct a kd-tree index:
    typedef KDTreeSingleIndexAdaptor<
      L2_Simple_Adaptor<double, PointCloud<double> > ,
      PointCloud<double>,
      3 /* dim */
      > my_kd_tree_t;

    my_kd_tree_t index(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf, to be tuned */) );

    index.buildIndex();

    std::vector<std::pair<long unsigned int,double> > ret_matches;

    nanoflann::SearchParams params;

    int maxNeigh = 0;
    for (int i = 0; i < *N; i++) {

      const double query_pt[3] = { cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z};

      const int nMatches = index.radiusSearch(&query_pt[0], *search_radius, ret_matches, params);
      maxNeigh = maxNeigh > nMatches ? maxNeigh : nMatches;

      cout << "radiusSearch(): radius=" << *search_radius << " -> " << nMatches << " matches\n";
      cout << "pt " << query_pt[0] << " " << query_pt[1] << " " << query_pt[2] << "\n";
      for (int j = 0; j < nMatches; j++) {
        cout << "idx["<< j << "]=" << ret_matches[j].first << " dist["<< j << "]=" << sqrt(ret_matches[j].second);
        cout << " " << cloud.pts[ret_matches[j].first].x << " " << cloud.pts[ret_matches[j].first].y << " " << cloud.pts[ret_matches[j].first].z;
        cout << " " << pow( cloud.pts[ret_matches[j].first].x - query_pt[0], 2) + pow( cloud.pts[ret_matches[j].first].y - query_pt[1], 2) + pow( cloud.pts[ret_matches[j].first].z - query_pt[2], 2);
        cout << endl;
      }
      cout << "\n";
    }
  }
}
