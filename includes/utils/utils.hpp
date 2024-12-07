#ifndef UTILS_HPP
#define UTILS_HPP

#include "../custom_cdt_class/custom_cdt_class.hpp"

typedef CGAL::Constrained_triangulation_plus_2<custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag>> CDT;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment;
typedef CDT::Point Point;
typedef CDT::Edge Edge;

namespace utils {
  
  struct FaceData {
    Point p1;
    Point p2;
    Point p3;

    FaceData(const Point& p1, const Point& p2, const Point& p3)
      : p1(p1), p2(p2), p3(p3) {}
  };

  class ant_parameters {
    public:
      double alpha;
      double beta;
      double xi;
      double psi;
      double lambda;
      double kappa;
      double L;

      ant_parameters(double alpha, double beta, double xi, double psi, double lambda, double kappa, double L);
  };

  class t_sp {
    public:
      double projection;
      double midpoint;
      double centroid;
      double circumcenter;
      double merge_obtuse;

      t_sp(double projection, double midpoint, double centroid, double circumcenter, double merge_obtuse);
  };

  class dt {
    public:
      double projection;
      double midpoint;
      double centroid;
      double circumcenter;
      double merge_obtuse;

      dt(double projection, double midpoint, double centroid, double circumcenter, double merge_obtuse);
  };

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
      std::list<CDT::Face_handle> affected_faces;

      obt_face(int count, CDT::Face_handle f);

      obt_face(int count, CDT::Face_handle f, std::list<CDT::Face_handle> affected_faces);
  };

  extern Polygon_2 region_boundary_polygon;
  
  // Check if a triangle is inside the region boundary
  bool is_triangle_inside_region_boundary(CDT::Face_handle f1);

  // Checks if a triangle has an obtuse angle
  bool has_obtuse_angle(CDT::Face_handle face);

  // Get the shared edge of two triangles
  Edge get_shared_edge(CDT &cdt, CDT::Face_handle f1, CDT::Face_handle neigh);

  // Count the number of obtuse triangles in the CDT
  int count_obtuse_triangles(CDT& cdt);

  // Find the obtuse vertex of a triangle
  int find_obtuse_vertex_id(CDT::Face_handle face);

  // Check if a polygon of 4 edges is convex
  bool is_convex(CDT& cdt, CDT::Face_handle f1, CDT::Face_handle f2);

  // Find the projection point of the obtuse vertex onto the opposite edge
  Point find_perpendicular_projection(CDT::Face_handle f, int obtuse_vertex_idx);

  // Check if a point is part of a constraint edge
  bool point_part_of_contrained_edge(CDT& cdt, Point p);

  // Check if a vertex is part of a constraint edge
  bool vertex_part_of_constrained_edge(CDT& cdt, CDT::Vertex_handle v);

  // Find an edge of a cdt by the points given
  void find_edge_by_points(CDT& cdt, Edge& edge, Point p1, Point p2);

  // Returns if two triangles are mergable
  bool are_mergable(CDT& cdt, CDT::Face_handle face, CDT::Face_handle neigh, Edge& shared_edge);

  // Check if a polygon is convex from the points given
  bool is_convex_polygon(const std::vector<Point>& points);

  // Remove the points given from the CDT
  void remove_points(CDT& cdt, std::set<CDT::Vertex_handle>& to_remove_points, std::vector<Point>& removed_points);

  // Insert a point at the center of a polygon from cdt
  Point calculate_centroid_coords(CDT& cdt, Point& p1, Point& p2, Point& p3,
                                    bool mergable_neigh1, 
                                    bool mergable_neigh2, 
                                    bool mergable_neigh3,
                                    CDT::Face_handle& neigh1, 
                                    CDT::Face_handle& neigh2, 
                                    CDT::Face_handle& neigh3);

  // Checks if the points are the same
  bool equal_points(Point a, Point b);

  // Check if two sets of points are the same
  bool equal_edges(Point a1, Point a2, Point b1, Point b2);

  // Get point from edge
  Point get_point_from_edge(Edge e, int point_number);

  // Get vertex from edge
  CDT::Vertex_handle get_vertex_from_edge(Edge e, int vertex_number);

  // Mark points to be removed
  void mark_points_to_remove(CDT& cdt,
                              Edge e, 
                              CDT::Face_handle neigh, 
                              std::vector<std::pair<Point, Point>>& edges_to_remove, 
                              std::set<CDT::Vertex_handle>& to_remove_points, 
                              CDT& polygon_cdt,
                              std::list<FaceData>& faces);

  // Checks if the faces are the same
  bool same_faces(FaceData f1, FaceData f2);

  // Checks if the face data is the same as the face
  bool same_FaceData_to_face(FaceData f1, CDT::Face_handle f2);

  // Get the face from face data
  CDT::Face_handle get_face_from_face_data(CDT& cdt, FaceData f);

  // Store the face data
  void store_face_data(CDT::Face_handle face, std::list<FaceData>& affected_faces);

  // Calculate the euclidean distance between two points
  double eucledean_distance(Point p1, Point p2);

  // Return the largest edge length of a triangle
  double largest_edge_length(CDT::Face_handle face);

  // Return the triagles height from the longest side
  double triangle_height_from_longest_side(CDT::Face_handle face);

  // Calculate the radius-to-height ratio
  double calculate_r_to_h(CDT::Face_handle face);

  // Check if there are 2 or more adjacent obtuse faces
  bool more_or_equal_to_2_adjacent_obtuse_faces(CDT& cdt, CDT::Face_handle face);

  // Calculate the posibility of a steiner method
  double calculate_posibility(double t, double h, double xi, double psi, double sum);

  // Update the dt value
  double update_dt(int obt_count, double dt, int steiner_counter, ant_parameters ant_params);

  // Update the pheromones
  void update_pheromones(t_sp& tsp, ant_parameters ant_params, dt Dt);
}

#endif // UTILS_HPP