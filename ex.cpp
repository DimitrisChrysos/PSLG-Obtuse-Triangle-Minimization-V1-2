#include "includes/custom_cdt_class/custom_cdt_class.hpp"
#include "includes/utils/utils.hpp"
#include "includes/read_write_file/read_write_file.hpp"
// #include <boost/json/src.hpp>
// #include <boost/json.hpp>
// #include <boost/property_tree/ptree.hpp>
// #include <boost/property_tree/json_parser.hpp>

typedef CGAL::Constrained_triangulation_plus_2<custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag>> CDT;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment;
typedef CDT::Point Point;
typedef CDT::Edge Edge;

using namespace utils;

// // Make flips in the CDT if possible and worth it
// void make_flips(CDT& cdt) {
//   for (const Edge& e : cdt.finite_edges()) {

//     // Get the edge points
//     CDT::Face_handle f1 = e.first;
//     int i = e.second;
//     CDT::Face_handle f2 = f1->neighbor(i);
//     Point a = get_point_from_edge(e, 1);
//     Point b = get_point_from_edge(e, 2);

//     // Test if the flip is possible or if it is worth doing
//     std::cout << "point a: " << a << " | point b: " << b << std::endl;
//     if (cdt.is_constrained(e) || !is_convex(cdt, f1, f2)) {
//       continue;
//     }


//     // If the flip is possible and worth it, do it
//     CDT copy(cdt);
//     if (test_the_flip(copy, a, b)) {
//       CGAL::draw(cdt);
//       cdt.tds().flip(f1, i);
//       cdt.remove_constrained_edge(f1, i);
//       CGAL::draw(cdt);
//     }
//   }
// }

// Insert the projection point of the obtuse vertex onto the opposite edge
obt_point insert_projection(CDT& cdt, CDT::Face_handle f1) {
  int obt_id = find_obtuse_vertex_id(f1);
  Point projection = find_perpendicular_projection(f1, obt_id);
  cdt.insert_no_flip(projection);

  obt_point ret(count_obtuse_triangles(cdt), projection);
  return ret;
}

// // Insert a point at the circumcenter of the triangle
// obt_point insert_circumcenter(CDT& cdt, CDT::Face_handle f1) {

//   // Calculate the circumcenter of the triangle
//   Point a = f1->vertex(0)->point();
//   Point b = f1->vertex(1)->point();
//   Point c = f1->vertex(2)->point();
//   Point pericenter = CGAL::circumcenter(a, b, c);

//   // Check if the point to be inserted will be inside the region_boundary_polygon
//   CDT::Face_handle located_face = cdt.locate(pericenter);
//   if (cdt.is_infinite(located_face) || !(CGAL::bounded_side_2(region_boundary_polygon.vertices_begin(), region_boundary_polygon.vertices_end(), pericenter) == CGAL::ON_BOUNDED_SIDE)) {    
//     obt_point ret(-1, pericenter);
//     return ret;
//   }

//   // CGAL::draw(cdt);

//   // Get the edge to be removed and
//   // check how many edges the main segment intersects
//   Point obtuse_point = f1->vertex(find_obtuse_vertex_id(f1))->point();
//   Segment main_segment(pericenter, obtuse_point);
//   Point intersect_point1;
//   Point intersect_point2;
//   CDT::Vertex_handle v1;
//   CDT::Vertex_handle v2;
//   Edge intersected_edge;
//   int count_intersect = 0;
//   for (const Edge& e : cdt.finite_edges()) {
//     Point edge_point1 = get_point_from_edge(e, 1);
//     Point edge_point2 = get_point_from_edge(e, 2);
//     Segment edge_segment(edge_point1, edge_point2);
//     if (CGAL::do_intersect(main_segment, edge_segment) && 
//         !equal_points(obtuse_point, edge_point1) && 
//         !equal_points(obtuse_point, edge_point2)) {

//       // Add an intersected edge
//       count_intersect++;
//       intersected_edge = e;
//       v1 = get_vertex_from_edge(e, 1);
//       v2 = get_vertex_from_edge(e, 2);
//       intersect_point1 = edge_point1;
//       intersect_point2 = edge_point2;

//       // If an interected edge is constrained
//       // or there are more than one intersecteed edges
//       // the method fails
//       if (cdt.is_constrained(e) || count_intersect > 1) {
//         obt_point ret(-1, pericenter);
//         return ret;
//       }
//     }
//   }

//   // Initialize the points to be removed
//   std::set<CDT::Vertex_handle> to_remove_points;
//   to_remove_points.insert(v1);
//   to_remove_points.insert(v2);

//   // Iterate over the faces of the points to be removed:
//   // Check if a point is part of a constrained edge
//   // If it is, remove the constrained edge and re-test if the point is part of a constrained edge
//   // If not, continue checking the next point
//   std::vector<std::pair<Point, Point>> false_removed_edges;
//   for (const CDT::Vertex_handle temp_vertex : to_remove_points) {

//     Point point = temp_vertex->point();

//     bool active = true;
//     while(active) {
//       Edge constrained_edge;
//       if (point_part_of_contrained_edge(cdt, point, false_removed_edges, constrained_edge)) {
//         CDT::Face_handle face = constrained_edge.first;
//         int index = constrained_edge.second;
//         cdt.remove_constrained_edge(face, index);
//       } else {
//         active = false;
//       }
//     }
//   }

//   // Remove the points
//   cdt.remove(v1);
//   cdt.remove(v2);
//   // CGAL::draw(cdt);

//   // Add the circumcenter point
//   cdt.insert_no_flip(pericenter);
//   // CGAL::draw(cdt);

//   // Add the removed edges as constrained edges
//   for (const auto& edge : false_removed_edges) {
//     cdt.insert_constraint(edge.first, edge.second);
//   }
//   CDT::Constraint_id cid = cdt.insert_constraint(intersect_point1, intersect_point2);
//   // CGAL::draw(cdt);

//   // Remove the constrained edge that needs to be removed
//   cdt.remove_constraint(cid);
//   // CGAL::draw(cdt);

//   // std::cout << "Finished insert_circumcenter\n";

//   obt_point ret(count_obtuse_triangles(cdt), pericenter);
//   return ret;
// }



// Insert a point at the circumcenter of the triangle
obt_point insert_circumcenter(CDT& cdt, CDT::Face_handle f1) {

  // Calculate the circumcenter of the triangle
  Point a = f1->vertex(0)->point();
  Point b = f1->vertex(1)->point();
  Point c = f1->vertex(2)->point();
  Point pericenter = CGAL::circumcenter(a, b, c);

  // Check if the point to be inserted will be inside the region_boundary_polygon
  CDT::Face_handle located_face = cdt.locate(pericenter);
  if (cdt.is_infinite(located_face) || !(CGAL::bounded_side_2(region_boundary_polygon.vertices_begin(), region_boundary_polygon.vertices_end(), pericenter) == CGAL::ON_BOUNDED_SIDE)) {    
    obt_point ret(-1, pericenter);
    return ret;
  }

  // Get the edge to be removed and
  // check how many edges the main segment intersects
  Point obtuse_point = f1->vertex(find_obtuse_vertex_id(f1))->point();
  Segment main_segment(pericenter, obtuse_point);
  Point intersect_point1;
  Point intersect_point2;
  CDT::Vertex_handle v1;
  CDT::Vertex_handle v2;
  Edge intersected_edge;
  int count_intersect = 0;
  for (const Edge& e : cdt.finite_edges()) {
    Point edge_point1 = get_point_from_edge(e, 1);
    Point edge_point2 = get_point_from_edge(e, 2);
    Segment edge_segment(edge_point1, edge_point2);
    if (CGAL::do_intersect(main_segment, edge_segment) && 
        !equal_points(obtuse_point, edge_point1) && 
        !equal_points(obtuse_point, edge_point2)) {

      // Add an intersected edge
      count_intersect++;
      intersected_edge = e;
      v1 = get_vertex_from_edge(e, 1);
      v2 = get_vertex_from_edge(e, 2);
      intersect_point1 = edge_point1;
      intersect_point2 = edge_point2;

      // If an interected edge is constrained
      // or there are more than one intersecteed edges
      // the method fails
      if (cdt.is_constrained(e) || count_intersect > 1) {
        obt_point ret(-1, pericenter);
        return ret;
      }
    }
  }

  // std::cout << "1..\n";

  // Remove the points able to be removed
  // If no points have the ability to be removed, return failure
  Edge constrained_edge;
  std::vector<std::pair<Point, Point>> removed_edges;
  bool fail = true;
  if (!point_part_of_contrained_edge(cdt, intersect_point1, removed_edges, constrained_edge)) {
    cdt.remove_no_flip(v1);
    fail = false;
  }
  if (!point_part_of_contrained_edge(cdt, intersect_point2, removed_edges, constrained_edge)) {
    cdt.remove_no_flip(v2);
    fail = false;
  }
  if (fail) {
    obt_point ret(-1, pericenter);
    return ret;
  }

  // std::cout << "2..\n";

  // Add the circumcenter point
  cdt.insert_no_flip(pericenter);

  // std::cout << "3..\n";


  // Add constrains to all the edges of the polygon, except the shared one
  std::vector<CDT::Constraint_id> cids;
  CDT::Face_handle neigh = f1->neighbor(intersected_edge.second);
  for (int i = 0; i < 3; ++i) {
    CDT::Vertex_handle start = neigh->vertex(CDT::ccw(i));
    CDT::Vertex_handle end = neigh->vertex(CDT::cw(i));
    Point start_point = start->point();
    Point end_point = end->point();
    if (!equal_edges(start_point, end_point, intersect_point1, intersect_point2)) {
      CDT::Constraint_id cid = cdt.insert_constraint(start_point, end_point);
      cids.push_back(cid);
    }
  }
  if (!equal_edges(a, b, intersect_point1, intersect_point2)) {
    CDT::Constraint_id cid = cdt.insert_constraint(a, b);
    cids.push_back(cid);
  }
  if (!equal_edges(a, c, intersect_point1, intersect_point2)) {
    CDT::Constraint_id cid = cdt.insert_constraint(a, c);
    cids.push_back(cid);
  }
  if (!equal_edges(b, c, intersect_point1, intersect_point2)) {
    CDT::Constraint_id cid = cdt.insert_constraint(b, c);
    cids.push_back(cid);
  }
  // CGAL::draw(cdt);

  // std::cout << "4..\n";

  // Remove the constrained edges
  for (const auto& cid : cids) {
    cdt.remove_constraint(cid);
  }
  // CGAL::draw(cdt);

  // std::cout << "5..\n";

  obt_point ret(count_obtuse_triangles(cdt), pericenter);
  return ret;
}


// Insert a point at the centroid of the triangle
obt_point insert_centroid(CDT& cdt, CDT::Face_handle f1) {

  // Calculate the centroid of the triangle
  Point a = f1->vertex(0)->point();
  Point b = f1->vertex(1)->point();
  Point c = f1->vertex(2)->point();
  Point centroid = CGAL::centroid(a, b, c);

  // Insert the centroid
  cdt.insert_no_flip(centroid);

  obt_point ret(count_obtuse_triangles(cdt), centroid);
  return ret;
}

// Merge triangles if possible
obt_face merge_obtuse(CDT& cdt, CDT::Face_handle f1) {

  obt_face ret(9999, f1);
  
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
  std::vector<std::pair<Point, Point>> edges_to_remove;
  std::vector<CDT::Face_handle> faces;
  faces.push_back(f1);
  if (mergable_neigh1) {
    mark_points_to_remove(cdt, e1, neigh1, edges_to_remove, to_remove_points, polygon_cdt, faces);
  }
  if (mergable_neigh2) {
    mark_points_to_remove(cdt, e2, neigh2, edges_to_remove, to_remove_points, polygon_cdt, faces);
  }
  if (mergable_neigh3) {
    mark_points_to_remove(cdt, e3, neigh3, edges_to_remove, to_remove_points, polygon_cdt, faces);
  }

  // If no points have the ability to be removed, return failure
  if (to_remove_points.size() == 0) {
    return ret;
  }

  // Check if the polygon is convex:
  // Iterate over the edges of the polygon_cdt
  // if there exist edges that are not constrained,
  // the polygon is not convex
  for (const Edge& e : polygon_cdt.finite_edges()) {
    if (!polygon_cdt.is_constrained(e)) {
      return ret;
    }
  }


  // Calculate centroid coords
  Point centroid = calculate_centroid_coords(cdt, p1, p2, p3, mergable_neigh1, mergable_neigh2, mergable_neigh3, neigh1, neigh2, neigh3);

  // Remove the points
  for (const CDT::Vertex_handle v : to_remove_points) {
    cdt.remove_no_flip(v);
  }



  // Add the centroid (or mean point) of the polygon
  cdt.insert_no_flip(centroid);
  cdt.insert_steiner_x_y(centroid.x(), centroid.y());

  // Add all the edges of the merged faces, except from the shared edges as constrains
  std::set<std::pair<Point, Point>> edges_made_constrained;
  for (const auto& face : faces) {
    for (int i = 0; i < 3; ++i) {
      Point temp_p1 = face->vertex((i + 1) % 3)->point();
      Point temp_p2 = face->vertex((i + 2) % 3)->point();

      bool to_remove = false;
      for (const auto& edge : edges_to_remove) {
        if (equal_edges(temp_p1, temp_p2, edge.first, edge.second)) {
          to_remove = true;
          break;
        }
      }
      if (!to_remove) {
        cdt.insert_constraint(temp_p1, temp_p2);
        edges_made_constrained.insert(std::make_pair(temp_p1, temp_p2));
        break;
      }
    }
  }


  // Remove the shared edges as constraints
  for (const auto& pair : edges_made_constrained) {
    Edge e;
    find_edge_by_points(cdt, e, pair.first, pair.second);
    cdt.remove_constrained_edge(e.first, e.second);
  }

  ret.obt_count = count_obtuse_triangles(cdt);
  return ret;
}

// Insert a point at the midpoint of the longest edge
obt_point insert_mid(CDT& cdt, CDT::Face_handle f1) {
  CDT copy(cdt);

  // Get the vertices of the triangle
  Point a = f1->vertex(0)->point();
  Point b = f1->vertex(1)->point();
  Point c = f1->vertex(2)->point();

  // Calculate the length of the edges
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

  cdt.insert_no_flip(midpoint);
  int obt = count_obtuse_triangles(cdt);
  obt_point ret = obt_point(obt, midpoint);
  return ret;
}

// // Choose the best method to insert a steiner point
// int steiner_insertion(CDT& cdt) {
//   int init_obtuse_count = count_obtuse_triangles(cdt);
//   Point a;
//   CDT::Face_handle f1;
//   obt_point best_steiner(9999, a);
//   obt_face of(9999, f1);
//   bool is_steiner_circumcenter = false;

//   // Iterate the faces of the cdt
//   for (CDT::Finite_faces_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++) {

//     if (!is_triangle_inside_region_boundary(f))
//       continue;

//     if (has_obtuse_angle(f)) {
      
//       CDT copy(cdt);
//       obt_point calc_insert_proj = insert_projection(copy, f);
//       if (best_steiner.obt_count >= calc_insert_proj.obt_count) {
//         best_steiner = calc_insert_proj;
//         is_steiner_circumcenter = false;
//       }

//       CDT copy1(cdt);
//       obt_point calc_insert_mid = insert_mid(copy1, f);
//       if (best_steiner.obt_count > calc_insert_mid.obt_count) {
//         is_steiner_circumcenter = false;
//         best_steiner = calc_insert_mid;
//       }

//       CDT copy2(cdt);
//       obt_point calc_insert_centr = insert_centroid(copy2, f);
//       if (best_steiner.obt_count > calc_insert_centr.obt_count) {
//         is_steiner_circumcenter = false;
//         best_steiner = calc_insert_centr;
//       }

//       CDT copy3(cdt);
//       obt_point calc_insert_circ = insert_circumcenter(copy3, f);
//       if (best_steiner.obt_count > calc_insert_circ.obt_count) {
//         is_steiner_circumcenter = true;
//         f1 = f;
//         best_steiner = calc_insert_circ;
//       }

//       CDT copy4(cdt);
//       obt_face temp = merge_obtuse(copy4, f);
//       if (temp.obt_count != -1 && temp.obt_count < of.obt_count) {
//         of = temp;
//       }
//     }
//   }
//   if (best_steiner.obt_count <= of.obt_count && best_steiner.obt_count <= count_obtuse_triangles(cdt)) {
//     if (is_steiner_circumcenter) {
//       insert_circumcenter(cdt, f1);
//     }
//     else {
//       cdt.insert_no_flip(best_steiner.insrt_pt);
//     }
//     cdt.insert_steiner_x_y(best_steiner.insrt_pt.x(), best_steiner.insrt_pt.y());
//     return 1;
//   }
//   else if (of.obt_count < best_steiner.obt_count && of.obt_count < count_obtuse_triangles(cdt)) {
//     merge_obtuse(cdt, of.face);
//     return 1;
//   }
//   return 0;
// }

enum class InsertionMethod {
  PROJECTION,
  MIDPOINT,
  CENTROID,
  CIRCUMCENTER,
  MERGE_OBTUSE,
  NONE
};

// Choose the best method to insert a steiner point
int steiner_insertion(CDT& cdt) {
  int init_obtuse_count = count_obtuse_triangles(cdt);
  Point starting_point;
  CDT::Face_handle starting_face;
  obt_point best_steiner(9999, starting_point);
  obt_face merge_face(9999, starting_face);
  CDT::Face_handle circumcenter_face;
  InsertionMethod best_method = InsertionMethod::NONE;

  // Iterate the faces of the cdt
  for (CDT::Finite_faces_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++) {

    if (!is_triangle_inside_region_boundary(f))
      continue;

    if (has_obtuse_angle(f)) {
      
      // std::cout << "1..\n";
      CDT copy4(cdt);
      obt_face temp_merge_face = merge_obtuse(copy4, f);
      if (best_steiner.obt_count > temp_merge_face.obt_count) {
        merge_face = temp_merge_face;
        best_method = InsertionMethod::MERGE_OBTUSE;
      }

      // std::cout << "2..\n";
      CDT copy3(cdt);
      obt_point calc_insert_circ = insert_circumcenter(copy3, f);
      if (calc_insert_circ.obt_count != -1 && best_steiner.obt_count > calc_insert_circ.obt_count) {
        circumcenter_face = f;
        best_steiner = calc_insert_circ;
        best_method = InsertionMethod::CIRCUMCENTER;
      }

      // std::cout << "3..\n";
      CDT copy2(cdt);
      obt_point calc_insert_centr = insert_centroid(copy2, f);
      if (best_steiner.obt_count > calc_insert_centr.obt_count) {
        best_steiner = calc_insert_centr;
        best_method = InsertionMethod::CENTROID;
      }

      // std::cout << "4..\n";
      CDT copy1(cdt);
      obt_point calc_insert_mid = insert_mid(copy1, f);
      if (best_steiner.obt_count > calc_insert_mid.obt_count) {
        best_steiner = calc_insert_mid;
        best_method = InsertionMethod::MIDPOINT;
      }

      // std::cout << "5..\n";
      CDT copy(cdt);
      obt_point calc_insert_proj = insert_projection(copy, f);
      if (best_steiner.obt_count >= calc_insert_proj.obt_count) {
        best_steiner = calc_insert_proj;
        best_method = InsertionMethod::PROJECTION;
      }
    }
  }
  if (best_method == InsertionMethod::PROJECTION || 
      best_method == InsertionMethod::MIDPOINT || 
      best_method == InsertionMethod::CENTROID) {
    cdt.insert_no_flip(best_steiner.insrt_pt);
    cdt.insert_steiner_x_y(best_steiner.insrt_pt.x(), best_steiner.insrt_pt.y());
    return 1;
  }
  else if (best_method == InsertionMethod::CIRCUMCENTER) {
    insert_circumcenter(cdt, circumcenter_face);
    cdt.insert_steiner_x_y(best_steiner.insrt_pt.x(), best_steiner.insrt_pt.y());
    return 1;
  }
  else if (best_method == InsertionMethod::MERGE_OBTUSE) {
    merge_obtuse(cdt, merge_face.face);
    return 1;
  }
  return 0;
}


// Create the region_boundary_polygon polygon
Polygon_2 make_region_boundary_polygon(std::list<int> region_boundary, std::vector<Point> points) {
  Polygon_2 region_boundary_polygon;
  for (int temp : region_boundary) {
    region_boundary_polygon.push_back(points[temp]);
  }
  return region_boundary_polygon;
}

using namespace read_write_file;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Wrong number of arguments\n";
    return 1;
  }
  
  namespace pt = boost::property_tree; // namespace alias
  pt::ptree root; // create a root node

  // Choose the file to read
  pt::read_json(argv[1], root); // read the json file
  
  // Read the json file
  std::string instance_uid = get_instance_uid(root);
  int num_points = get_num_points(root);
  std::list<int> points_x = get_points_x(root);
  std::list<int> points_y = get_points_y(root);
  std::list<int> region_boundary = get_region_boundary(root);
  std::string num_constraints = get_num_constraints(root);
  std::list<std::pair<int, int>> additional_constraints = get_additional_constraints(root, region_boundary);
  // std::string method = get_method(root);
  // std::list<std::pair<std::string, double>> parameters = get_parameters(root);
  // bool delaunay = get_delaunay(root);

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

  // Create the region boundary polygon
  region_boundary_polygon = make_region_boundary_polygon(region_boundary, points);

  // Count the obtuse triangles
  // std::cout << "Before flips | obt_triangles: " << count_obtuse_triangles(cdt) << std::endl;
  CGAL::draw(cdt);

  // Make flips
  // make_flips(cdt);
  // CGAL::draw(cdt);
  // std::cout << "After flips | obt_triangles: " << count_obtuse_triangles(cdt) << std::endl;

  // Insert Steiner points
  int steps = 0;
  int consec_insertions = 0;
  int count_obt = 0;
  for (int i = 0 ; i < 200 ; i++) {
    if (count_obtuse_triangles(cdt) == 0) {
      break;
    }
    steps += steiner_insertion(cdt);
    if (count_obt == count_obtuse_triangles(cdt)) {
      consec_insertions++;
    }
    else {
      consec_insertions = 0;
    }
    count_obt = count_obtuse_triangles(cdt);
    // if (consec_insertions > 40) break;
    std::cout << "After try " << i << " to insert Steiner | obt_triangles: " << count_obtuse_triangles(cdt) << std::endl;
  }
  std::cout << "After " << steps << " Steiner Insertions | obt_triangles: " << count_obtuse_triangles(cdt) << std::endl;
  CGAL::draw(cdt);

  write_output(cdt, points);
  
  return 0;
}