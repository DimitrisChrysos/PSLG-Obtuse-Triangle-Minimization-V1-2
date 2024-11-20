#ifndef UTILS_HPP
#define UTILS_HPP

#include "../custom_cdt_class/custom_cdt_class.hpp"

typedef CGAL::Constrained_triangulation_plus_2<custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag>> CDT;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment;
typedef CDT::Point Point;
typedef CDT::Edge Edge;

namespace utils {
  class obt_point {
    public:
      int obt_count;
      Point insrt_pt;

      obt_point(int count, Point pt);
  };

  class obt_face {
    public:
      int obt_count;
      CDT::Face_handle face;

      obt_face(int count, CDT::Face_handle f);
  };

  extern Polygon_2 region_boundary_polygon;
  
  // Check if a triangle is inside the region boundary
  bool is_triangle_inside_region_boundary(CDT::Face_handle f1);

  // Checks if a triangle has an obtuse angle
  bool has_obtuse_angle(CDT::Face_handle face);

  // Get the shared edge of two triangles
  Edge get_shared_edge(CDT &cdt, CDT::Face_handle f1, CDT::Face_handle neigh);

  // Count the number of obtuse triangles in the CDT
  int count_obtuse_triangles(CDT cdt);

  // Find the obtuse vertex of a triangle
  int find_obtuse_vertex_id(CDT::Face_handle face);

  // Check if a polygon of 4 edges is convex
  bool is_convex(CDT& cdt, CDT::Face_handle f1, CDT::Face_handle f2);

  // Find the projection point of the obtuse vertex onto the opposite edge
  Point find_perpendicular_projection(CDT::Face_handle f, int obtuse_vertex_idx);

  // Check if a point is part of a constraint edge
  bool point_part_of_contrained_edge(CDT& cdt, Point p);

  // Returns if two triangles are mergable
  bool are_mergable(CDT& cdt, CDT::Face_handle face, CDT::Face_handle neigh, Edge& shared_edge);

  // Check if a polygon is convex from the points given
  bool is_convex_polygon(const std::vector<Point>& points);

  // Test if the flip is possible
  bool test_the_flip(CDT& cdt, Point v1, Point v2);

  // Remove the points given from the CDT
  void remove_points(CDT& cdt, std::set<CDT::Vertex_handle>& to_remove_points, std::vector<Point>& removed_points);

  // Checks if the points are the same
  bool equal_points(Point a, Point b);
}

#endif // UTILS_HPP