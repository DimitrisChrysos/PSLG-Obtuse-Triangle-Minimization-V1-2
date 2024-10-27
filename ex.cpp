#include "includes/utils/utils.hpp"

// Make flips in the CDT if possible and worth it
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

// Insert the projection point of the obtuse vertex onto the opposite edge
obt_point insert_projection(CDT& cdt, CDT::Face_handle f1) {
  int obt_id = find_obtuse_vertex_id(f1);
  Point projection = find_perpendicular_projection(f1, obt_id);
  cdt.insert_no_flip(projection);

  obt_point ret(count_obtuse_triangles(cdt), projection);
  return ret;
}

// Insert a point at the circumcenter of the triangle
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

  // Check if the polygon is convex:
  // Iterate over the edges of the polygon_cdt
  // if there exist edges that are not constrained,
  // the polygon is not convex
  for (const Edge& e : polygon_cdt.finite_edges()) {
    if (!polygon_cdt.is_constrained(e)) {
      return ret;
    }
  }

  // Remove the points
  std::vector<Point> removed_points;
  remove_points(cdt, to_remove_points, removed_points);

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
  cdt.insert_no_flip(centroid);
  cdt.insert_steiner_x_y(centroid.x(), centroid.y());
  
  // Add the removed edge as a constrained edge
  std::vector<CDT::Constraint_id> constraint_ids;
  for (const auto& edge : edges_made_constrained) {
    CDT::Constraint_id cid = cdt.insert_constraint(edge.first, edge.second);
    constraint_ids.push_back(cid);
  }

  // Remove the constraint edge
  for (const auto& cid : constraint_ids) {
    cdt.remove_constraint(cid);
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
  Point mid;

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

// Choose the best method to insert a steiner point
int steiner_insertion(CDT& cdt) {
  int init_obtuse_count = count_obtuse_triangles(cdt);
  Point a;
  CDT::Face_handle f1;
  obt_point best_steiner(9999, a);
  obt_face of(9999, f1);

  // Iterate the faces of the cdt
  for (CDT::Finite_faces_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++) {

    if (!is_triangle_inside_region_boundary(f))
      continue;

    if (has_obtuse_angle(f)) {
      
      CDT copy(cdt);
      obt_point calc_insert_proj = insert_projection(copy, f);
      if (best_steiner.obt_count >= calc_insert_proj.obt_count) {
        best_steiner = calc_insert_proj;
      }

      CDT copy1(cdt);
      obt_point calc_insert_mid = insert_mid(copy1, f);
      if (best_steiner.obt_count > calc_insert_mid.obt_count) {
        best_steiner = calc_insert_mid;
      }

      CDT copy2(cdt);
      obt_point calc_insert_centr = insert_centroid(copy2, f);
      if (best_steiner.obt_count > calc_insert_centr.obt_count) {
        best_steiner = calc_insert_centr;
      }

      CDT copy3(cdt);
      obt_point calc_insert_circ = insert_circumcenter(copy3, f);
      if (best_steiner.obt_count > calc_insert_circ.obt_count) {
        best_steiner = calc_insert_circ;
      }

      CDT copy4(cdt);
      obt_face temp = merge_obtuse(copy4, f);
      if (temp.obt_count != -1 && temp.obt_count < of.obt_count) {
        of = temp;
      }
    }
  }
  if (best_steiner.obt_count <= of.obt_count && best_steiner.obt_count <= count_obtuse_triangles(cdt)) {
    cdt.insert_no_flip(best_steiner.insrt_pt);
    cdt.insert_steiner_x_y(best_steiner.insrt_pt.x(), best_steiner.insrt_pt.y());
    return 1;
  }
  else if (of.obt_count < best_steiner.obt_count && of.obt_count < count_obtuse_triangles(cdt)) {
    merge_obtuse(cdt, of.face);
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

int main() {
  
  namespace pt = boost::property_tree; // namespace alias
  pt::ptree root; // create a root node

  // Choose the file to read
  pt::read_json("input.json", root); // read the json file
  // pt::read_json("test_instances/instance_test_12.json", root); // read the json file
  // pt::read_json("tests/instance_7.json", root); // read the json file
  
  // Read the json file
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

  // Create the region boundary polygon
  region_boundary_polygon = make_region_boundary_polygon(region_boundary, points);

  // Count the obtuse triangles
  std::cout << "Before flips | obt_triangles: " << count_obtuse_triangles(cdt) << std::endl;
  CGAL::draw(cdt);

  // Make flips
  make_flips(cdt);
  std::cout << "After flips | obt_triangles: " << count_obtuse_triangles(cdt) << std::endl;

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
    if (consec_insertions > 40) break;
    std::cout << "After try to insert Steiner | obt_triangles: " << count_obtuse_triangles(cdt) << std::endl;
  }
  std::cout << "After " << steps << " Steiner Insertions | obt_triangles: " << count_obtuse_triangles(cdt) << std::endl;
  CGAL::draw(cdt);

  write_output(cdt, points);
  
  return 0;
}