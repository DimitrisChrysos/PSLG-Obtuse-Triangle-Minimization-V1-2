#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2.h>
// #include <CGAL/Constraint_id.h>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS // optional in native ubuntu, removes a warning in wsl
#ifndef CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H

#define CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H



#include <CGAL/Constrained_Delaunay_triangulation_2.h>



template <class Gt, class Tds = CGAL::Default, class Itag = CGAL::Default>

class Custom_Constrained_Delaunay_triangulation_2

    : public CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag> {

public:

    using Base = CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>;

    using typename Base::Face_handle;

    using typename Base::Point;

    using typename Base::Vertex_handle;

    using typename Base::Locate_type;



    // Constructors

    Custom_Constrained_Delaunay_triangulation_2(const Gt& gt = Gt())

        : Base(gt) {}



    Custom_Constrained_Delaunay_triangulation_2(typename Base::List_constraints& lc, const Gt& gt = Gt())

        : Base(lc, gt) {}



    template <class InputIterator>

    Custom_Constrained_Delaunay_triangulation_2(InputIterator it, InputIterator last, const Gt& gt = Gt())

        : Base(it, last, gt) {}



    // New insert method without flips

    Vertex_handle insert_no_flip(const Point& a, Face_handle start = Face_handle()) {

        // Call Ctr::insert without flip_around

        Vertex_handle va = this->Base::Ctr::insert(a, start); // Directly call Ctr::insert from the base

        return va;

    }



    // Another insert method with known location

    Vertex_handle insert_no_flip(const Point& a, Locate_type lt, Face_handle loc, int li) {

        Vertex_handle va = this->Base::Ctr::insert(a, lt, loc, li); // Directly call Ctr::insert from the base

        return va;

    }

    bool my_is_flippable(const typename Base::Edge& e) {
      Face_handle f1 = e.first; 
      int i = e.second; 
      Face_handle f2 = f1->neighbor(i);
      if (this->is_infinite(f1) || this->is_infinite(f2)) {
        return false;
      }
      else if (this->is_constrained(e)) {
        return false;
      }

      return true;
    }
};

#endif // CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H


typedef CGAL:: Exact_predicates_exact_constructions_kernel K;
typedef CGAL:: Exact_predicates_tag Itag;
// typedef Custom_Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<Custom_Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag>> CDT;
typedef CGAL::Polygon_2<K> Polygon_2;

// typedef K::Point_2 Point;
typedef K::Segment_2 Segment;

typedef CDT::Point Point;
typedef CDT::Edge Edge;

//
class obt_point {
  public:
    int obt_count;
    Point insrt_pt;

    obt_point(int count, Point pt) {
      obt_count = count;
      insrt_pt = pt;
    }
};

class obt_face {
  public:
    int obt_count;
    CDT::Face_handle face;

    obt_face(int count, CDT::Face_handle f) {
      obt_count = count;
      face = f;
    }
};
//

bool is_triangle_inside_region_boundary(Polygon_2& region_boundary_polygon, CDT::Face_handle f1);
Polygon_2 region_boundary_polygon;

bool has_obtuse_angle(CDT::Face_handle face) {
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

int count_obtuse_triangles(CDT cdt) {
  int count = 0;
  for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); fit++) {
    CDT::Face_handle face = fit;
    
    // Check if the face is inside the region boundary
    if (!is_triangle_inside_region_boundary(region_boundary_polygon, face))
      continue;

    if (has_obtuse_angle(face)) {
      count++;
    }
  }
  return count;
}

int find_obtuse_vertex_id(CDT::Face_handle face) {
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

Edge get_shared_edge(CDT &cdt, CDT::Face_handle f1, CDT::Face_handle neigh);

// check if a polygon of 4 edges is convex
bool is_convex(CDT& cdt, CDT::Face_handle f1, CDT::Face_handle f2) {

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

  // // Print them
  // for (auto p : shared_points) {
  //   std::cout << "shared_points: " << p << std::endl;
  // }
  // for (auto p : not_shared_points) {
  //   std::cout << "not_shared_points: " << p << std::endl;
  // }

  std::vector<Point> polygon_points = {}; 
  // Take one shared point and them remove it
  Point shared_point = *shared_points.begin();
  shared_points.erase(shared_points.begin());
  polygon_points.push_back(shared_point);
  // std::cout << "a: " << shared_point << std::endl;
  // Take the first not shared point and remove it
  Point not_shared_point = *not_shared_points.begin();
  not_shared_points.erase(not_shared_points.begin());
  polygon_points.push_back(not_shared_point);
  // std::cout << "b: " << not_shared_point << std::endl;
  // Take the second shared point and remove it
  shared_point = *shared_points.begin();
  shared_points.erase(shared_points.begin());
  polygon_points.push_back(shared_point);
  // std::cout << "c: " << shared_point << std::endl;
  // Take the second not shared point and remove it
  not_shared_point = *not_shared_points.begin();
  not_shared_points.erase(not_shared_points.begin());
  polygon_points.push_back(not_shared_point);
  // std::cout << "d: " << not_shared_point << std::endl;

  // // print the polygon points
  // for (auto p : polygon_points) {
  //   std::cout << "polygon_points: " << p << std::endl;
  // }

  return CGAL::is_convex_2(polygon_points.begin(), polygon_points.end());
}

bool test_the_flip(CDT& cdt, Point v1, Point v2) {
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

void make_flips(CDT& cdt) {
  int count = 0;
  for (const Edge& e : cdt.finite_edges()) {

    // Get the edge points
    CDT::Face_handle f1 = e.first;
    int i = e.second;
    CDT::Face_handle f2 = f1->neighbor(i);
    auto v1 = f1->vertex((i+1)%3);
    auto v2 = f1->vertex((i+2)%3);


    // Test is the flip possible or if it is worth doing
    bool do_flip = test_the_flip(cdt, v1->point(), v2->point());

    // If the flip is possible and worth it, do it
    if (do_flip) {
      cdt.tds().flip(f1, i);
      count++;
    }
  }
}

Point find_perpendicular_projection(CDT::Face_handle f, int obtuse_vertex_idx) {
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

obt_point insert_projection(CDT& cdt, CDT::Face_handle f1) {
  int obt_id = find_obtuse_vertex_id(f1);
  Point projection = find_perpendicular_projection(f1, obt_id);
  cdt.insert_no_flip(projection);

  obt_point ret(count_obtuse_triangles(cdt), projection);
  return ret;
}

Point compute_incenter(CDT::Face_handle f) {
  // Get the vertices of the face
  Point A = f->vertex(0)->point();
  Point B = f->vertex(1)->point();
  Point C = f->vertex(2)->point();

  // Calculate the lengths of the sides
  K::FT a = CGAL::squared_distance(B, C); // Side opposite to A
  K::FT b = CGAL::squared_distance(A, C); // Side opposite to B
  K::FT c = CGAL::squared_distance(A, B); // Side opposite to C

  // Calculate the incenter using the weighted average formula
  K::FT sum = a + b + c;
  Point incenter = CGAL::barycenter(A, a, B, b, C, c);

  return incenter;
}

obt_point insert_incenter(CDT& cdt, CDT::Face_handle f1) {
  Point incenter = compute_incenter(f1);
  cdt.insert_no_flip(incenter);

  obt_point ret(count_obtuse_triangles(cdt), incenter);
  return ret;
}

obt_point insert_circumcenter(CDT& cdt, CDT::Face_handle f1) {

  // Calculate the circumcenter of the triangle
  Point a = f1->vertex(0)->point();
  Point b = f1->vertex(1)->point();
  Point c = f1->vertex(2)->point();
  Point pericenter = CGAL::circumcenter(a, b, c);

  // Check if the inserted vertex is inside the convex hull
  CDT::Face_handle located_face = cdt.locate(pericenter);
  if (!cdt.is_infinite(located_face)) {
    cdt.insert_no_flip(pericenter);
  }

  // return count_obtuse_triangles(copy);
  obt_point ret(count_obtuse_triangles(cdt), pericenter);
  return ret;
}

obt_point insert_centroid(CDT& cdt, CDT::Face_handle f1) {

  // Calculate the centroid of the triangle
  Point a = f1->vertex(0)->point();
  Point b = f1->vertex(1)->point();
  Point c = f1->vertex(2)->point();
  Point centroid = CGAL::centroid(a, b, c);

  // Insert the centroid
  cdt.insert_no_flip(centroid);


  // count_obtuse_triangles(copy);
  obt_point ret(count_obtuse_triangles(cdt), centroid);
  return ret;
}



Edge get_shared_edge(CDT &cdt, CDT::Face_handle f1, CDT::Face_handle neigh) {
  // Get the index of the shared edge
  int edge_index = f1->index(neigh);

  // Get the vertices of the shared edge
  Point p1 = f1->vertex((edge_index + 1) % 3)->point();
  Point p2 = f1->vertex((edge_index + 2) % 3)->point();
  // std::cout << "1. Edge: " << p1 << " - " << p2 << std::endl;

  for (const Edge& e : cdt.finite_edges()) {
    CDT::Face_handle face = e.first;
    int index = e.second;
    Point edgeP1 = face->vertex((index + 1) % 3)->point();
    Point edgeP2 = face->vertex((index + 2) % 3)->point();

    if ((edgeP1 == p1 && edgeP2 == p2) || (edgeP1 == p2 && edgeP2 == p1)) {
      return e;
    }
  }
}

// Check if a point is part of a constraint edge
bool point_part_of_contrained_edge(CDT& cdt, Point p) {
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

// If two triangles are mergable
bool are_mergable(CDT& cdt, CDT::Face_handle face, CDT::Face_handle neigh, Edge& shared_edge) {
  
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

void remove_points(CDT& cdt, std::set<CDT::Vertex_handle>& to_remove_points, std::vector<Point>& removed_points) {
  
  // Remove the points
  for (CDT::Vertex_handle v : to_remove_points) {
    removed_points.push_back(v->point()); 
    cdt.remove(v);
  }
}

bool is_convex_polygon(const std::vector<Point>& points) {
  return CGAL::is_convex_2(points.begin(), points.end());
}


obt_face merge_obtuse(CDT& cdt, CDT::Face_handle f1) {

  // CDT copy(cdt);
  obt_face ret(-1, f1);

  // Get the vertices of the triangle
  Point p1 = f1->vertex(0)->point();
  Point p2 = f1->vertex(1)->point();
  Point p3 = f1->vertex(2)->point();

  // Get the neighbors of the triangle
  CDT::Face_handle neigh1 = f1->neighbor(0);
  CDT::Face_handle neigh2 = f1->neighbor(1);
  CDT::Face_handle neigh3 = f1->neighbor(2);

  // Get the shared edges of the triangle
  Edge e1 = get_shared_edge(cdt, f1, neigh1);
  Edge e2 = get_shared_edge(cdt, f1, neigh2);
  Edge e3 = get_shared_edge(cdt, f1, neigh3);

  // To check if polygon convex later
  CDT polygon_cdt;
  for (int i = 0; i < 3; ++i) {
    Point temp_p1 = f1->vertex((i + 1) % 3)->point();
    Point temp_p2 = f1->vertex((i + 2) % 3)->point();
    polygon_cdt.insert_constraint(temp_p1, temp_p2);
  }


  // Check if any of the neighbors are mergable
  bool mergable_neigh1 = are_mergable(cdt, f1, neigh1, e1);
  bool mergable_neigh2 = are_mergable(cdt, f1, neigh2, e2);
  bool mergable_neigh3 = are_mergable(cdt, f1, neigh3, e3);  

  // If none of the neighbors are mergable, return
  if (!mergable_neigh1 && !mergable_neigh2 && !mergable_neigh3) {
    return ret;
  }

  // If a neighbor is mergable, prepare the points to be removed
  std::set<CDT::Vertex_handle> to_remove_points;
  std::vector<std::pair<Point, Point>> edges_made_constrained;
  std::vector<CDT::Face_handle> faces;
  if (mergable_neigh1) {
    // To check convexity later
    for (int i = 0; i < 3; ++i) {
      Point temp_p1 = neigh1->vertex((i + 1) % 3)->point();
      Point temp_p2 = neigh1->vertex((i + 2) % 3)->point();
      polygon_cdt.insert_constraint(temp_p1, temp_p2);
    }

    // Mark points to be removed
    CDT::Vertex_handle v1 = e1.first->vertex((e1.second+1)%3);
    CDT::Vertex_handle v2 = e1.first->vertex((e1.second+2)%3);
    Point a = v1->point();
    Point b = v2->point();
    edges_made_constrained.push_back(std::make_pair(a, b));
    if(!point_part_of_contrained_edge(cdt, a)) to_remove_points.insert(v1);
    if(!point_part_of_contrained_edge(cdt, b)) to_remove_points.insert(v2);
    faces.push_back(neigh1);
  }
  if (mergable_neigh2) {
    // To check convexity later
    for (int i = 0; i < 3; ++i) {
      Point temp_p1 = neigh2->vertex((i + 1) % 3)->point();
      Point temp_p2 = neigh2->vertex((i + 2) % 3)->point();
      polygon_cdt.insert_constraint(temp_p1, temp_p2);
    }

    // Mark points to be removed
    CDT::Vertex_handle v1 = e2.first->vertex((e2.second+1)%3);
    CDT::Vertex_handle v2 = e2.first->vertex((e2.second+2)%3);
    Point a = v1->point();
    Point b = v2->point();
    edges_made_constrained.push_back(std::make_pair(a, b));
    if(!point_part_of_contrained_edge(cdt, a)) to_remove_points.insert(v1);
    if(!point_part_of_contrained_edge(cdt, b)) to_remove_points.insert(v2);
    faces.push_back(neigh2);
  }
  if (mergable_neigh3) {
    // To check convexity later
    for (int i = 0; i < 3; ++i) {
      Point temp_p1 = neigh3->vertex((i + 1) % 3)->point();
      Point temp_p2 = neigh3->vertex((i + 2) % 3)->point();
      polygon_cdt.insert_constraint(temp_p1, temp_p2);
    }

    // Mark points to be removed
    CDT::Vertex_handle v1 = e3.first->vertex((e3.second+1)%3);
    CDT::Vertex_handle v2 = e3.first->vertex((e3.second+2)%3);
    Point a = v1->point();
    Point b = v2->point();
    edges_made_constrained.push_back(std::make_pair(a, b));
    if(!point_part_of_contrained_edge(cdt, a)) to_remove_points.insert(v1);
    if(!point_part_of_contrained_edge(cdt, b)) to_remove_points.insert(v2);
    faces.push_back(neigh3);
  }


  // If to_remove_points is empty, return
  if (to_remove_points.empty()) {
    return ret;
  }


  // CGAL::draw(polygon_cdt);


  // Check if the polygon is convex:
  // Iterate over the edges of the polygon_cdt
  // if there exist edges that are not constrained,
  // the polygon is not convex
  for (const Edge& e : polygon_cdt.finite_edges()) {
    if (!polygon_cdt.is_constrained(e)) {
      return ret;
    }
  }

  // CGAL::draw(polygon_cdt);

  // Remove the points
  std::vector<Point> removed_points;
  remove_points(cdt, to_remove_points, removed_points);


  // CGAL::draw(cdt);


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
  std::cout << "centroid: " << centroid << std::endl;
  cdt.insert_no_flip(centroid);
  
  // CGAL::draw(cdt);

  // Add the removed edge as a constrained edge
  std::vector<CDT::Constraint_id> constraint_ids;
  for (const auto& edge : edges_made_constrained) {
    CDT::Constraint_id cid = cdt.insert_constraint(edge.first, edge.second);
    constraint_ids.push_back(cid);
  }

  // CGAL::draw(cdt);


  // Remove the constraint edge
  for (const auto& cid : constraint_ids) {
    cdt.remove_constraint(cid);
  }

  // CGAL::draw(cdt);



  ret.obt_count = count_obtuse_triangles(cdt);
  return ret;
}


obt_point insert_mid(CDT& cdt, CDT::Face_handle f1) {
  CDT copy(cdt);

  // Calculate the circumcenter of the triangle
  Point a = f1->vertex(0)->point();
  Point b = f1->vertex(1)->point();
  Point c = f1->vertex(2)->point();
  Point mid;

  K::FT l0 = CGAL::squared_distance(a, b);
  K::FT l1 = CGAL::squared_distance(a, c);
  K::FT l2 = CGAL::squared_distance(b, c);

  // Calculate the midpoint of edge with the longest length
  Point midpoint;
  if (l0 >= l1 && l0 >= l2) {
    midpoint = CGAL::midpoint(a, b);
  }
  else if (l1 >= l2) {
    midpoint = CGAL::midpoint(a, c);
  }
  else {
    midpoint = CGAL::midpoint(b, c);
  }

  // std::cout << "midpoint: " << midpoint << " | points: a:" << a << " b: " << b << std::endl;
  cdt.insert_no_flip(midpoint);
  int obt = count_obtuse_triangles(cdt);
  obt_point ret = obt_point(obt, midpoint);
  return ret;
}

bool is_triangle_inside_region_boundary(Polygon_2& region_boundary_polygon, CDT::Face_handle f1) {

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

void steiner_insertion(CDT& cdt, Polygon_2& region_boundary_polygon) {
  int init_obtuse_count = count_obtuse_triangles(cdt);
  std::cout << "Initial obtuse count: " << init_obtuse_count << std::endl;
  Point a;
  CDT::Face_handle f1;
  obt_point best_steiner(9999, a);
  obt_face of(9999, f1);

  // Iterate the faces of the cdt
  for (CDT::Finite_faces_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++) {

    if (!is_triangle_inside_region_boundary(region_boundary_polygon, f))
      continue;

    if (has_obtuse_angle(f)) {
      
      // CDT copy(cdt);
      // // Insert the circumcenter if possible
      // obt_point calc_insert_mid = insert_mid(copy, f);
      // if (best_steiner.obt_count > calc_insert_mid.obt_count) {
      //   best_steiner = calc_insert_mid;
      // }
      // // std::cout << "Obtuse triangles after inserting the edge midpoint: " << count_obtuse_triangles(copy) << std::endl;

      // CDT copy1(cdt);
      // // Insert the circumcenter if possible
      // obt_point calc_insert_centr = insert_centroid(copy1, f);
      // if (best_steiner.obt_count > calc_insert_centr.obt_count) {
      //   best_steiner = calc_insert_centr;
      // }
      // // std::cout << "Obtuse triangles after inserting the centroid: " << count_obtuse_triangles(copy1) << std::endl;

      // CDT copy2(cdt);
      // // Insert the circumcenter if possible
      // obt_point calc_insert_circ = insert_circumcenter(copy2, f);
      // if (best_steiner.obt_count > calc_insert_circ.obt_count) {
      //   best_steiner = calc_insert_circ;
      // }
      // // std::cout << "Obtuse triangles after inserting the circumcenter: " << count_obtuse_triangles(copy2) << std::endl;


      // CDT copy3(cdt);
      // // Insert the circumcenter if possible
      // obt_point calc_insert_inc = insert_incenter(copy3, f);
      // if (best_steiner.obt_count > calc_insert_inc.obt_count) {
      //   best_steiner = calc_insert_inc;
      // }

      CDT copy4(cdt);
      // Insert the circumcenter if possible
      obt_point calc_insert_proj = insert_projection(copy4, f);
      if (best_steiner.obt_count > calc_insert_proj.obt_count) {
        best_steiner = calc_insert_proj;
      }

      CDT copy5(cdt);
      obt_face temp = merge_obtuse(copy5, f);
      if (temp.obt_count != -1 && temp.obt_count < of.obt_count) {
        of = temp;
      }
    }
  }
  if (best_steiner.obt_count <= of.obt_count && best_steiner.obt_count < count_obtuse_triangles(cdt)) {
    cdt.insert_no_flip(best_steiner.insrt_pt);
    // na kanw kai insert to steiner_x_y
  }
  else if (of.obt_count < best_steiner.obt_count && of.obt_count < count_obtuse_triangles(cdt)) {
    merge_obtuse(cdt, of.face);
    std::cout << "mpaino sto merge gia to kanoniko polygono!\n";
  }

  std::cout << "Final obtuse count: " << count_obtuse_triangles(cdt) << std::endl;
}


// Read JSON functions
std::string get_instance_uid(boost::property_tree::ptree root) {
  return root.get<std::string>("instance_uid");
}

int get_num_points(boost::property_tree::ptree root) {
  return root.get<int>("num_points", 0);;
}

std::list<int> get_points_x(boost::property_tree::ptree root) {
  std::list<int> points_x;
  for (boost::property_tree::ptree::value_type &point_x : root.get_child("points_x")) {
    points_x.push_back(point_x.second.get_value<int>());
  }
  return points_x;
}

std::list<int> get_points_y(boost::property_tree::ptree root) {
  std::list<int> points_y;
  for (boost::property_tree::ptree::value_type &point_y : root.get_child("points_y")) {
    points_y.push_back(point_y.second.get_value<int>());
  }
  return points_y;
}

std::list<int> get_region_boundary(boost::property_tree::ptree root) {
  std::list<int> region_boundary;
  for (boost::property_tree::ptree::value_type &temp : root.get_child("region_boundary")) {
    region_boundary.push_back(temp.second.get_value<int>());
  }
  return region_boundary;
}

std::string get_num_constraints(boost::property_tree::ptree root) {
  return root.get<std::string>("num_constraints");
}

std::list<std::pair<int, int>> get_additional_constraints(boost::property_tree::ptree root, std::list<int> region_boundary) {
  std::list<std::pair<int, int>> additional_constraints;
  for (boost::property_tree::ptree::value_type &row : root.get_child("additional_constraints")) {
    auto it = row.second.begin();
    int first = it->second.get_value<int>();
    ++it;
    int second = it->second.get_value<int>();
    additional_constraints.push_back(std::make_pair(first, second));
  }

  // Add region_boundary to additional_constraints
  int prev = -1;
  int first;
  for (const auto &temp : region_boundary) {
    if (prev == -1) {
      prev = temp;
      first = temp;
      continue;
    }
    additional_constraints.push_back(std::make_pair(prev, temp));
    prev = temp;
  }
  additional_constraints.push_back(std::make_pair(prev, first));

  return additional_constraints;
}

Polygon_2 make_region_boundary_polygon(std::list<int> region_boundary, std::vector<Point> points) {
  
  // Create region_boundary_polygon polygon
  Polygon_2 region_boundary_polygon;
  for (int temp : region_boundary) {
    std::cout << "temp: " << temp << std::endl;
    region_boundary_polygon.push_back(points[temp]);
  }

  return region_boundary_polygon;
}

int main() {
  
  // Read the json file
  namespace pt = boost::property_tree; // namespace alias
  pt::ptree root; // create a root node
  // pt::read_json("input.json", root); // read the json file
  pt::read_json("test_instances/instance_test_14.json", root); // read the json file
  std::string instance_uid = get_instance_uid(root);
  int num_points = get_num_points(root);
  std::list<int> points_x = get_points_x(root);
  std::list<int> points_y = get_points_y(root);
  std::list<int> region_boundary = get_region_boundary(root);
  std::string num_constraints = get_num_constraints(root);
  std::list<std::pair<int, int>> additional_constraints = get_additional_constraints(root, region_boundary);


  // Create the Constrained Delaunay Triangulation (CDT)
  CDT cdt;

  // Define the points from the PSLG (x, y coordinates) and insert them into the CDT
  std::vector<Point> points = {};
  auto it_y = points_y.begin();
  for (const auto& p : points_x) {
    int point_x = p;
    int point_y = *it_y;
    points.push_back(Point(point_x, point_y));
    cdt.insert(Point(point_x, point_y));
    it_y++;
  }

  // Add the constrained edges from additional_constraints
  for (const auto &constraint : additional_constraints) {
    cdt.insert_constraint(points[constraint.first], points[constraint.second]);
  }

  region_boundary_polygon = make_region_boundary_polygon(region_boundary, points);
  for (auto it = region_boundary_polygon.vertices_begin(); it != region_boundary_polygon.vertices_end(); ++it) {
    std::cout << "(" << it->x() << ", " << it->y() << ")" << std::endl;
  }

  CGAL::draw(cdt);


  // Count the obtuse triangles
  std::cout << "Number of obtuse triangles before the flips: " << count_obtuse_triangles(cdt) << std::endl;

  // Make flips
  make_flips(cdt);

  CGAL::draw(cdt);

  // Insert Steiner points
  std::cout << "\n\n\nSteiner points insertion:\n";
  for (int i = 0 ; i < 1 ; i++) {
    steiner_insertion(cdt, region_boundary_polygon);
  }

  // Count the obtuse triangles
  std::cout << "Number of obtuse triangles after the flips: " << count_obtuse_triangles(cdt) << std::endl;

  // Draw the triangulation using CGAL's draw function
  CGAL::draw(cdt);

  return 0;
}