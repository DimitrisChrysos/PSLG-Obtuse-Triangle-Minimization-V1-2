#include "utils.hpp"
#include <random>

namespace utils {
  Polygon_2 region_boundary_polygon;
}

utils::ant_parameters::ant_parameters(double alpha, double beta, double xi, double psi, double lambda, double kappa, double L) {
  this->alpha = alpha;
  this->beta = beta;
  this->xi = xi;
  this->psi = psi;
  this->lambda = lambda;
  this->kappa = kappa;
  this->L = L;
}

utils::t_sp::t_sp(double projection, double midpoint, double centroid, double circumcenter, double merge_obtuse) {
  this->projection = projection;
  this->midpoint = midpoint;
  this->centroid = centroid;
  this->circumcenter = circumcenter;
  this->merge_obtuse = merge_obtuse;
}

utils::dt::dt(double projection, double midpoint, double centroid, double circumcenter, double merge_obtuse) {
  this->projection = projection;
  this->midpoint = midpoint;
  this->centroid = centroid;
  this->circumcenter = circumcenter;
  this->merge_obtuse = merge_obtuse;
}

utils::obt_point::obt_point(int count, Point pt) {
  this->obt_count = count;
  this->insrt_pt = pt;
}

utils::obt_face::obt_face(int count, CDT::Face_handle f) {
  this->obt_count = count;
  this->face = f;
}

utils::obt_face::obt_face(int count, CDT::Face_handle f, std::list<CDT::Face_handle> affected_faces) 
  : obt_count(count), face(f), affected_faces(std::move(affected_faces)) {}

// Check if a triangle is inside the region boundary
bool utils::is_triangle_inside_region_boundary(CDT::Face_handle f1) {
  // Get the vertices of the triangle
  Point p1 = f1->vertex(0)->point();
  Point p2 = f1->vertex(1)->point();
  Point p3 = f1->vertex(2)->point();

  // Get the centroid of the triangle
  Point centroid = CGAL::centroid(p1, p2, p3);
  if (CGAL::bounded_side_2(region_boundary_polygon.vertices_begin(), region_boundary_polygon.vertices_end(), centroid) == CGAL::ON_BOUNDED_SIDE) {
    return true;
  }
  return false;
}

// Checks if a triangle has an obtuse angle
bool utils::has_obtuse_angle(CDT::Face_handle face) {
  // Get the vertices of the triangle
  Point p1 = face->vertex(0)->point();
  Point p2 = face->vertex(1)->point();
  Point p3 = face->vertex(2)->point();

  // Check if any angle of the triangle is obtuse
  if (CGAL::angle(p1,p2,p3) == CGAL::OBTUSE ||
      CGAL::angle(p2,p3,p1) == CGAL::OBTUSE ||
      CGAL::angle(p3,p1,p2) == CGAL::OBTUSE) {
    return true;
  }
  return false;
}

// Get the shared edge of two triangles
Edge utils::get_shared_edge(CDT &cdt, CDT::Face_handle f1, CDT::Face_handle neigh) {

  // Get the index of the shared edge
  int edge_index = f1->index(neigh);

  // Get the vertices of the shared edge
  Point p1 = f1->vertex((edge_index + 1) % 3)->point();
  Point p2 = f1->vertex((edge_index + 2) % 3)->point();

  const Edge ed;
  for (const Edge& e : cdt.finite_edges()) {
    CDT::Face_handle face = e.first;
    int index = e.second;
    Point edgeP1 = face->vertex((index + 1) % 3)->point();
    Point edgeP2 = face->vertex((index + 2) % 3)->point();

    if ((edgeP1 == p1 && edgeP2 == p2) || (edgeP1 == p2 && edgeP2 == p1)) {
      return e;
    }
  }

  return ed;
}

// Count the number of obtuse triangles in the CDT
int utils::count_obtuse_triangles(CDT& cdt) {
  int count = 0;
  for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); fit++) {
    CDT::Face_handle face = fit;
    
    // Check if the face is inside the region boundary
    if (!is_triangle_inside_region_boundary(face))
      continue;

    if (has_obtuse_angle(face)) {
      count++;
    }
  }
  return count;
}

// Find the obtuse vertex of a triangle
int utils::find_obtuse_vertex_id(CDT::Face_handle face) {
  Point p1 = face->vertex(0)->point();
  Point p2 = face->vertex(1)->point();
  Point p3 = face->vertex(2)->point();

  // Check if any angle of the triangle is obtuse
  if (CGAL::angle(p1,p2,p3) == CGAL::OBTUSE) {
    return 1;
  }
  else if (CGAL::angle(p2,p3,p1) == CGAL::OBTUSE) {
    return 2;
  }
  else if (CGAL::angle(p3,p1,p2) == CGAL::OBTUSE) {
    return 0;
  }
  return -1;
}

// Check if a polygon of 4 edges is convex
bool utils::is_convex(CDT& cdt, CDT::Face_handle f1, CDT::Face_handle f2) {

  // Get the points of f1
  Point p1 = f1->vertex(0)->point();
  Point p2 = f1->vertex(1)->point();
  Point p3 = f1->vertex(2)->point();

  // Get the points of f2
  Point p4 = f2->vertex(0)->point();
  Point p5 = f2->vertex(1)->point();
  Point p6 = f2->vertex(2)->point();

  // Not Shared Points -> not_shared_points
  std::set<Point> not_shared_points = {p1, p2, p3};
  auto it = not_shared_points.find(p4);
  if (it != not_shared_points.end()) {
      not_shared_points.erase(it);
  }
  else {
    not_shared_points.insert(p4);
  }
  it = not_shared_points.find(p5);
  if (it != not_shared_points.end()) {
      not_shared_points.erase(it);
  }
  else {
    not_shared_points.insert(p5);
  }
  it = not_shared_points.find(p6);
  if (it != not_shared_points.end()) {
      not_shared_points.erase(it);
  }
  else {
    not_shared_points.insert(p6);
  }

  // Shared Points -> "shared_points"
  std::set<Point> temp_points = {p1, p2, p3};
  std::set<Point> shared_points = {};
  it = temp_points.find(p4);
  if (it != temp_points.end()) {
    shared_points.insert(p4);
  }
  it = temp_points.find(p5);
  if (it != temp_points.end()) {
    shared_points.insert(p5);
  }
  it = temp_points.find(p6);
  if (it != temp_points.end()) {
    shared_points.insert(p6);
  }

  // Create a polygon with the points of the two triangles
  std::vector<Point> polygon_points = {}; 
  // Take one shared point and them remove it
  Point shared_point = *shared_points.begin();
  shared_points.erase(shared_points.begin());
  polygon_points.push_back(shared_point);
  // Take the first not shared point and remove it
  Point not_shared_point = *not_shared_points.begin();
  not_shared_points.erase(not_shared_points.begin());
  polygon_points.push_back(not_shared_point);
  // Take the second shared point and remove it
  shared_point = *shared_points.begin();
  shared_points.erase(shared_points.begin());
  polygon_points.push_back(shared_point);
  // Take the second not shared point and remove it
  not_shared_point = *not_shared_points.begin();
  not_shared_points.erase(not_shared_points.begin());
  polygon_points.push_back(not_shared_point);

  // Return if the polygon is convex
  return CGAL::is_convex_2(polygon_points.begin(), polygon_points.end());
}

// Find the projection point of the obtuse vertex onto the opposite edge
Point utils::find_perpendicular_projection(CDT::Face_handle f, int obtuse_vertex_idx) {
  // Get the vertices of the face
  Point A = f->vertex(0)->point();
  Point B = f->vertex(1)->point();
  Point C = f->vertex(2)->point();
  
  Point obtuse_vertex, vertex1, vertex2;
  
  // Identify which vertex has the obtuse angle
  if (obtuse_vertex_idx == 0) {
      obtuse_vertex = A;
      vertex1 = B;
      vertex2 = C;
  } else if (obtuse_vertex_idx == 1) {
      obtuse_vertex = B;
      vertex1 = A;
      vertex2 = C;
  } else {
      obtuse_vertex = C;
      vertex1 = A;
      vertex2 = B;
  }

  // Project the obtuse vertex onto the opposite edge (vertex1, vertex2)
  Segment opposite_edge(vertex1, vertex2);

  // Project the obtuse vertex onto the opposite edge
  Point projection = opposite_edge.supporting_line().projection(obtuse_vertex);

  return projection;
}

// Check if a point is part of a constraint edge
bool utils::point_part_of_contrained_edge(CDT& cdt, Point p) {
  for (const Edge& e : cdt.finite_edges()) {
    if (cdt.is_constrained(e)) {
      Point edge_point1 = get_point_from_edge(e, 1);
      Point edge_point2 = get_point_from_edge(e, 2);

      if ((equal_points(edge_point1, p) || equal_points(edge_point2, p))) {
        return true;
      }
    }
  }
  return false;
}

// Check if a vertex is part of a constraint edge
bool utils::vertex_part_of_constrained_edge(CDT& cdt, CDT::Vertex_handle v) {
  int counter=0;
  for (auto it = cdt.incident_edges(v); it != nullptr; ++it) {
    if (cdt.is_constrained(*it)) {
      return true; // Vertex is part of a constrained edge
    }
    counter++;
    if (counter > 100) {
      return true;
    }
  }
  return false;
}

// Find an edge of a cdt by the points given
void utils::find_edge_by_points(CDT& cdt, Edge& edge, Point p1, Point p2) {
  for (const Edge& e : cdt.finite_edges()) {
    Point edge_point1 = get_point_from_edge(e, 1);
    Point edge_point2 = get_point_from_edge(e, 2);
    if ((equal_points(edge_point1, p1) && equal_points(edge_point2, p2)) || 
        (equal_points(edge_point1, p2) && equal_points(edge_point2, p1))) {
      edge = e;
      return;
    }
  }  
}

// Returns if two triangles are mergable
bool utils::are_mergable(CDT& cdt, CDT::Face_handle face, CDT::Face_handle neigh, Edge& shared_edge) {
  
  // If the neighbor is not obtused or their shared edge is constrained return false
  if (!has_obtuse_angle(neigh) || cdt.is_constrained(shared_edge))
    return false;

  // If both points of the shared_edge are part of a constrained edge return false
  CDT::Vertex_handle v1 = get_vertex_from_edge(shared_edge, 1);
  CDT::Vertex_handle v2 = get_vertex_from_edge(shared_edge, 2);
  if (vertex_part_of_constrained_edge(cdt, v1) || vertex_part_of_constrained_edge(cdt, v2)) {
    return false;
  }
  Point p1 = get_point_from_edge(shared_edge, 1);
  Point p2 = get_point_from_edge(shared_edge, 2);
  if (point_part_of_contrained_edge(cdt, p1) ||
      point_part_of_contrained_edge(cdt, p2)) {
    return false;
  }


  return true;
}

// Check if a polygon is convex from the points given
bool utils::is_convex_polygon(const std::vector<Point>& points) {
  return CGAL::is_convex_2(points.begin(), points.end());
}

// Insert a point at the center of a polygon from cdt
Point utils::calculate_centroid_coords(CDT& cdt, Point& p1, Point& p2, Point& p3,
                                    bool mergable_neigh1, 
                                    bool mergable_neigh2, 
                                    bool mergable_neigh3,
                                    CDT::Face_handle& neigh1, 
                                    CDT::Face_handle& neigh2, 
                                    CDT::Face_handle& neigh3) {

  // Add the centroid (or mean point) of the polygon
  std::vector<Point> points = {
    p1, p2, p3
  };
  if (mergable_neigh1) {
    points.push_back(neigh1->vertex(0)->point());
    points.push_back(neigh1->vertex(1)->point());
    points.push_back(neigh1->vertex(2)->point());
  }
  if (mergable_neigh2) {
    points.push_back(neigh2->vertex(0)->point());
    points.push_back(neigh2->vertex(1)->point());
    points.push_back(neigh2->vertex(2)->point());
  }
  if (mergable_neigh3) {
    points.push_back(neigh3->vertex(0)->point());
    points.push_back(neigh3->vertex(1)->point());
    points.push_back(neigh3->vertex(2)->point());
  }

  K::FT mean_x = 0;
  K::FT mean_y = 0;
  for (const auto& point : points) {
    mean_x += point.x();
    mean_y += point.y();
  }
  K::FT points_size = points.size();
  mean_x /= points_size;
  mean_y /= points_size;
  Point centroid(mean_x, mean_y);
  return centroid;
}

// Checks if the points are the same
bool utils::equal_points(Point a, Point b) {
  if (a.x() == b.x() && a.y() == b.y()) {
    return true;
  }
  return false;
}

// Check if two sets of points are the same
bool utils::equal_edges(Point a1, Point a2, Point b1, Point b2) {
  if ((equal_points(a1, b1) && equal_points(a2, b2)) || 
      (equal_points(a1, b2) && equal_points(a2, b1))) {
    return true;
  }
  return false;
}

// Get point from edge
Point utils::get_point_from_edge(Edge e, int point_number) {
  CDT::Face_handle face = e.first;
  int index = e.second;
  return face->vertex((index + point_number) % 3)->point();
}

// Get vertex from edge
CDT::Vertex_handle utils::get_vertex_from_edge(Edge e, int vertex_number) {
  CDT::Face_handle face = e.first;
  int index = e.second;
  return face->vertex((index + vertex_number) % 3);
}

void utils::mark_points_to_remove(CDT& cdt, 
                                  Edge e, 
                                  CDT::Face_handle neigh, 
                                  std::vector<std::pair<Point, Point>>& edges_to_remove, 
                                  std::set<CDT::Vertex_handle>& to_remove_points, 
                                  CDT& polygon_cdt,
                                  std::list<FaceData>& faces) {

  CDT::Vertex_handle v1 = get_vertex_from_edge(e, 1);
  CDT::Vertex_handle v2 = get_vertex_from_edge(e, 2);
  Point a = v1->point();
  Point b = v2->point();
  Edge constrained_edge;
  to_remove_points.insert(v1);
  to_remove_points.insert(v2);
  edges_to_remove.push_back(std::make_pair(a, b));
  for (int i = 0; i < 3; ++i) {
    Point temp_p1 = neigh->vertex((i + 1) % 3)->point();
    Point temp_p2 = neigh->vertex((i + 2) % 3)->point();
    polygon_cdt.insert_constraint(temp_p1, temp_p2);
  }
  FaceData f = {neigh->vertex(0)->point(), neigh->vertex(1)->point(), neigh->vertex(2)->point()};
  faces.push_back(f);
}

// Checks if the faces are the same
bool utils::same_faces(utils::FaceData f1, utils::FaceData f2) {

  // If a combination of the points are the same return true
  if (equal_points(f1.p1, f2.p1) && equal_points(f1.p2, f2.p2) && equal_points(f1.p3, f2.p3)) {
    return true;
  }
  if (equal_points(f1.p1, f2.p1) && equal_points(f1.p2, f2.p3) && equal_points(f1.p3, f2.p2)) {
    return true;
  }
  if (equal_points(f1.p1, f2.p2) && equal_points(f1.p2, f2.p1) && equal_points(f1.p3, f2.p3)) {
    return true;
  }
  if (equal_points(f1.p1, f2.p2) && equal_points(f1.p2, f2.p3) && equal_points(f1.p3, f2.p1)) {
    return true;
  }
  if (equal_points(f1.p1, f2.p3) && equal_points(f1.p2, f2.p1) && equal_points(f1.p3, f2.p2)) {
    return true;
  }
  if (equal_points(f1.p1, f2.p3) && equal_points(f1.p2, f2.p2) && equal_points(f1.p3, f2.p1)) {
    return true;
  }
  return false;
}

bool utils::same_FaceData_to_face(FaceData f1, CDT::Face_handle f2) {
  Point p1 = f2->vertex(0)->point();
  Point p2 = f2->vertex(1)->point();
  Point p3 = f2->vertex(2)->point();

  // If a combination of the points are the same return true
  if (equal_points(f1.p1, p1) && equal_points(f1.p2, p2) && equal_points(f1.p3, p3)) {
    return true;
  }
  if (equal_points(f1.p1, p1) && equal_points(f1.p2, p3) && equal_points(f1.p3, p2)) {
    return true;
  }
  if (equal_points(f1.p1, p2) && equal_points(f1.p2, p1) && equal_points(f1.p3, p3)) {
    return true;
  }
  if (equal_points(f1.p1, p2) && equal_points(f1.p2, p3) && equal_points(f1.p3, p1)) {
    return true;
  }
  if (equal_points(f1.p1, p3) && equal_points(f1.p2, p1) && equal_points(f1.p3, p2)) {
    return true;
  }
  if (equal_points(f1.p1, p3) && equal_points(f1.p2, p2) && equal_points(f1.p3, p1)) {
    return true;
  }
  return false;
}

// Get the face from face data
CDT::Face_handle utils::get_face_from_face_data(CDT& cdt, FaceData f) {
  for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); fit++) {
    CDT::Face_handle face = fit;
    if (same_FaceData_to_face(f, face)) {
      return face;
    }
  }
  return cdt.finite_faces_begin();
}

// Store the face data
void utils::store_face_data(CDT::Face_handle face, std::list<FaceData>& affected_faces) {
  Point p1 = face->vertex(0)->point();
  Point p2 = face->vertex(1)->point();
  Point p3 = face->vertex(2)->point();
  affected_faces.emplace_back(p1, p2, p3);
}

// Calculate the euclidean distance between two points
double utils::eucledean_distance(Point p1, Point p2) {
  double p1x = CGAL::to_double(p1.x());
  double p1y = CGAL::to_double(p1.y());
  double p2x = CGAL::to_double(p2.x());
  double p2y = CGAL::to_double(p2.y());
  return std::sqrt(std::pow(p1x - p2x, 2) + std::pow(p1y - p2y, 2));
}

// Return the largest edge length of a triangle
double utils::largest_edge_length(CDT::Face_handle face) {
  // Get the vertices of the triangle
  Point a = face->vertex(0)->point();
  Point b = face->vertex(1)->point();
  Point c = face->vertex(2)->point();

  // Calculate the length of the edges
  double l0 = eucledean_distance(a, b);
  double l1 = eucledean_distance(a, c);
  double l2 = eucledean_distance(b, c);

  // Return largest edge length 
  if (l0 >= l1 && l0 >= l2) {
    return l0;
  }
  else if (l1 >= l2) {
    return l1;
  }
  else {
    return l2;
  }
}

// Return the triagles height from the longest side
double utils::triangle_height_from_longest_side(CDT::Face_handle face) {
  Point p1 = face->vertex(0)->point();
  Point p2 = face->vertex(1)->point();
  Point p3 = face->vertex(2)->point();
  if (CGAL::angle(p1,p2,p3) == CGAL::OBTUSE) {
    int obt_id = 1;
    Point projection = find_perpendicular_projection(face, obt_id);
    return eucledean_distance(p2, projection);
  }
  else if (CGAL::angle(p2,p3,p1) == CGAL::OBTUSE) {
    int obt_id = 2;
    Point projection = find_perpendicular_projection(face, obt_id);
    return eucledean_distance(p3, projection);
  }
  else if (CGAL::angle(p3,p1,p2) == CGAL::OBTUSE) {
    int obt_id = 0;
    Point projection = find_perpendicular_projection(face, obt_id);
    return eucledean_distance(p1, projection);
  }
  return -1;
}

// Calculate the radius-to-height ratio
double utils::calculate_r_to_h(CDT::Face_handle face) {
  return (double)(( (double)largest_edge_length(face) / (double)2) / (double)triangle_height_from_longest_side(face));
}

// Check if there are 2 or more adjacent obtuse faces
bool utils::more_or_equal_to_2_adjacent_obtuse_faces(CDT& cdt, CDT::Face_handle face) {
  int obtuse_count = 0;
  for (int i = 0; i < 3; i++) {
    CDT::Face_handle neigh = face->neighbor(i);
    Edge e = get_shared_edge(cdt, face, neigh);

    if (has_obtuse_angle(neigh) && are_mergable(cdt, face, neigh, e)) {
      obtuse_count++;
    }
  }
  return obtuse_count >= 2;
}

// Calculate the posibility of a steiner method
double utils::calculate_posibility(double t, double h, double xi, double psi, double sum) {
  return (double)( std::pow(t, xi) * std::pow(h, psi) ) / sum;
}

// Update the dt value
double utils::update_dt(int obt_count, double dt, int steiner_counter, ant_parameters ant_params) {
  dt = dt + ((double)1 / ( (double)1 + ant_params.alpha*(double)obt_count + ant_params.beta*(double)steiner_counter ));
  return dt;
}

// Update the pheromones
void utils::update_pheromones(t_sp& tsp, ant_parameters ant_params, dt Dt) {
  tsp.projection = ((double)1 - ant_params.lambda)*tsp.projection + Dt.projection;
  tsp.midpoint = ((double)1 - ant_params.lambda)*tsp.midpoint + Dt.midpoint;
  tsp.centroid = ((double)1 - ant_params.lambda)*tsp.centroid + Dt.centroid;
  tsp.circumcenter = ((double)1 - ant_params.lambda)*tsp.circumcenter + Dt.circumcenter;
  tsp.merge_obtuse = ((double)1 - ant_params.lambda)*tsp.merge_obtuse + Dt.merge_obtuse;
}