#include "utils.hpp"

namespace utils {
  Polygon_2 region_boundary_polygon;
}

utils::obt_point::obt_point(int count, Point pt) {
  obt_count = count;
  insrt_pt = pt;
}

utils::obt_face::obt_face(int count, CDT::Face_handle f) {
  obt_count = count;
  face = f;
}

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
int utils::count_obtuse_triangles(CDT cdt) {
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
      CDT::Face_handle face = e.first;
      int index = e.second;
      Point edge_point1 = face->vertex((index + 1) % 3)->point();
      Point edge_point2 = face->vertex((index + 2) % 3)->point();

      if ((edge_point1.x() == p.x() && edge_point1.y() == p.y()) 
        || (edge_point2.x() == p.x() && edge_point2.y() == p.y()))
        return true;
    }
  }
  return false;
}

// Returns if two triangles are mergable
bool utils::are_mergable(CDT& cdt, CDT::Face_handle face, CDT::Face_handle neigh, Edge& shared_edge) {
  
  // If the neighbor is not obtused or their shared edge is constrained return false
  if (!has_obtuse_angle(neigh) || cdt.is_constrained(shared_edge))
    return false;

  // Get the points of the shared_edge
  CDT::Face_handle face1 = shared_edge.first;
  int index = shared_edge.second;
  Point edge_point1 = face1->vertex((index + 1) % 3)->point();
  Point edge_point2 = face1->vertex((index + 2) % 3)->point();

  // If both of the points of the shared_edge are part of a constrained edge the triangles are not mergable
  if (point_part_of_contrained_edge(cdt, edge_point1) && point_part_of_contrained_edge(cdt, edge_point2))
    return false;

  return true;
}

// Check if a polygon is convex from the points given
bool utils::is_convex_polygon(const std::vector<Point>& points) {
  return CGAL::is_convex_2(points.begin(), points.end());
}

// Test if the flip is possible
bool utils::test_the_flip(CDT& cdt, Point v1, Point v2) {
  CDT copy(cdt);
  for (const Edge& e : copy.finite_edges()) {

    // Get the edge points
    CDT::Face_handle f1 = e.first; 
    int i = e.second; 
    CDT::Face_handle f2 = f1->neighbor(i);
    auto p1 = f1->vertex((i+1)%3);
    auto p2 = f1->vertex((i+2)%3);
    Point p1_point = p1->point();
    Point p2_point = p2->point();

    // Check if the edge is the one we are looking for
    if (!(v1.x() == p1_point.x() && v1.y() == p1_point.y() && v2.x() == p2_point.x() && v2.y() == p2_point.y())) {
      continue;
    }

    // Check if the polygon created is convex
    if (!is_convex(copy, f1, f2)) {
      return false;
    };

    // Check if the edge is flippable
    if (!copy.my_is_flippable(e)) {
      return false;
    }

    // Check if the triangles formed by the edge have obtuse angles
    int obt = 0;
    if (has_obtuse_angle(f1)) obt++;
    if (has_obtuse_angle(f2)) obt++;

    // If triangles have obtuse angles, make the flip
    if (obt) {

      // Make the flip
      copy.tds().flip(f1, i);

      // Check if the triangles formed by the edge have obtuse angles after the flip
      int obt2 = 0;
      if (has_obtuse_angle(f1)) obt2++;
      if (has_obtuse_angle(f2)) obt2++;

      // If the number of obtuse angles after the flip are less, return true
      if (obt2 < obt) {
        return true;
      }
    }

    return false;
  }

  return false;
}

// Remove the points given from the CDT
void utils::remove_points(CDT& cdt, std::set<CDT::Vertex_handle>& to_remove_points, std::vector<Point>& removed_points) {
  
  // Remove the points
  for (CDT::Vertex_handle v : to_remove_points) {
    std::cout << "Removing point: " << v->point() << std::endl;
    removed_points.push_back(v->point());
    cdt.remove(v);
  }
}

// Checks if the points are the same
bool utils::equal_points(Point a, Point b) {
  if (a.x() == b.x() && a.y() == b.y()) {
    return true;
  }
  return false;
}