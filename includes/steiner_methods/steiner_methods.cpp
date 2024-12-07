#include "steiner_methods.hpp"

steiner_methods::InsertionMethod steiner_methods::choose_random_steiner_method() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, static_cast<int>(InsertionMethod::MERGE_OBTUSE));
  int selected = dis(gen);
  if (selected == static_cast<int>(InsertionMethod::PROJECTION)) {
    return InsertionMethod::PROJECTION;
  }
  else if (selected == static_cast<int>(InsertionMethod::MIDPOINT)) {
    return InsertionMethod::MIDPOINT;
  }
  else if (selected == static_cast<int>(InsertionMethod::CENTROID)) {
    return InsertionMethod::CENTROID;
  }
  else if (selected == static_cast<int>(InsertionMethod::CIRCUMCENTER)) {
    return InsertionMethod::CIRCUMCENTER;
  }
  else if (selected == static_cast<int>(InsertionMethod::MERGE_OBTUSE)) {
    return InsertionMethod::MERGE_OBTUSE;
  }

  return InsertionMethod::NONE;
}

// Insert the projection point of the obtuse vertex onto the opposite edge
obt_point steiner_methods::insert_projection(CDT& cdt, CDT::Face_handle f1) {
  int obt_id = find_obtuse_vertex_id(f1);
  Point projection = find_perpendicular_projection(f1, obt_id);
  cdt.insert_no_flip(projection);
  obt_point ret(count_obtuse_triangles(cdt), projection);
  return ret;
}

// Insert a point at the centroid of the triangle
obt_point steiner_methods::insert_centroid(CDT& cdt, CDT::Face_handle f1) {

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

// Insert a point at the midpoint of the longest edge
obt_point steiner_methods::insert_mid(CDT& cdt, CDT::Face_handle f1) {
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

// Insert a point at the circumcenter of the triangle
obt_face steiner_methods::insert_circumcenter(CDT& cdt, FaceData toReplaceFace) {

  // Find the face from FaceData
  CDT::Face_handle f1 = get_face_from_face_data(cdt, toReplaceFace);

  // Initialize the return value
  obt_face ret(9999, f1);

  // Calculate the circumcenter of the triangle
  Point a = f1->vertex(0)->point();
  Point b = f1->vertex(1)->point();
  Point c = f1->vertex(2)->point();

  Edge bc = std::make_pair(f1, 0);
  Edge ac = std::make_pair(f1, 1);
  Edge ab = std::make_pair(f1, 2);
  
  // Get the circumcenter of the triangle
  Point pericenter = CGAL::circumcenter(a, b, c);

  // Check if the point to be inserted will be inside the region_boundary_polygon
  CDT::Face_handle located_face = cdt.locate(pericenter);
  if (cdt.is_infinite(located_face) || !(CGAL::bounded_side_2(region_boundary_polygon.vertices_begin(), region_boundary_polygon.vertices_end(), pericenter) == CGAL::ON_BOUNDED_SIDE)) {    
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
    v1 = get_vertex_from_edge(e, 1);
    v2 = get_vertex_from_edge(e, 2);
    Point edge_point1 = v1->point();
    Point edge_point2 = v2->point();
    Segment edge_segment(edge_point1, edge_point2);
    if (CGAL::do_intersect(main_segment, edge_segment) && 
        !equal_points(obtuse_point, edge_point1) && 
        !equal_points(obtuse_point, edge_point2)) {

      // Add an intersected edge
      count_intersect++;
      intersected_edge = e;
      intersect_point1 = v1->point();
      intersect_point2 = v2->point();
    }
  }

  // If the interected edge is constrained
  // or there are more than one intersected edges
  // the method fails
  if (cdt.is_constrained(intersected_edge) || count_intersect > 1) {
    return ret;
  }

  // Get a copy of the cdt
  CDT copy(cdt);
  CDT::Face_handle f1_copy = get_face_from_face_data(copy, toReplaceFace);

  // Get the neighbor of the face
  CDT::Face_handle neigh_of_f1 = f1->neighbor(intersected_edge.second);
  Point neigh_a = neigh_of_f1->vertex(0)->point();
  Point neigh_b = neigh_of_f1->vertex(1)->point();
  Point neigh_c = neigh_of_f1->vertex(2)->point();

  // Remove the points if possible, else fail
  bool fail = true;
  if (point_part_of_contrained_edge(cdt, a) ||
      point_part_of_contrained_edge(cdt, b) ||
      point_part_of_contrained_edge(cdt, c) ||
      point_part_of_contrained_edge(cdt, neigh_a) ||
      point_part_of_contrained_edge(cdt, neigh_b) ||
      point_part_of_contrained_edge(cdt, neigh_c)) {
    return ret;
  }
  if (vertex_part_of_constrained_edge(cdt, v1) || vertex_part_of_constrained_edge(cdt, v2)) {
    return ret;
  }
  cdt.remove_no_flip(v1);
  cdt.remove_no_flip(v2);

  // Add the circumcenter point
  cdt.insert_no_flip(pericenter);
  cdt.insert_steiner_x_y(pericenter.x(), pericenter.y());


  // Add constrains to all the edges of the polygon created, except the shared one
  std::set<Edge> edges_to_remove;
  std::vector<CDT::Constraint_id> cids;
  CDT::Face_handle neigh = f1_copy->neighbor(intersected_edge.second);
  for (int i = 0; i < 3; ++i) {
    Point start_point = neigh->vertex((i + 1) % 3)->point();
    Point end_point = neigh->vertex((i + 2) % 3)->point();
    Edge e = std::make_pair(neigh, i);

    if (!equal_edges(start_point, end_point, intersect_point1, intersect_point2)) {
      if (cdt.is_constrained(e)) { // If the edge is constrained, skip
        continue;
      }
      cdt.insert_constraint(start_point, end_point);
      edges_to_remove.insert(e);
    }
  }
  if (!equal_edges(a, b, intersect_point1, intersect_point2)){//} && !cdt.is_constrained(ab)) {
    cdt.insert_constraint(a, b);
    edges_to_remove.insert(ab);
  }
  if (!equal_edges(a, c, intersect_point1, intersect_point2)){// && !cdt.is_constrained(ac)) {
    cdt.insert_constraint(a, c);
    edges_to_remove.insert(ac);
  }
  if (!equal_edges(b, c, intersect_point1, intersect_point2)){// && !cdt.is_constrained(bc)) {
    cdt.insert_constraint(b, c);
    edges_to_remove.insert(bc);
  }

  // Remove the constrained edges
  for (const auto& edge : edges_to_remove) {
    cdt.remove_constrained_edge(edge.first, edge.second);
  }


  ret.obt_count = count_obtuse_triangles(cdt);
  std::list<CDT::Face_handle> affected_faces;
  affected_faces.push_back(f1);
  affected_faces.push_back(neigh);
  ret.affected_faces = affected_faces;
  return ret;
}

// Merge triangles if possible
obt_face steiner_methods::merge_obtuse(CDT& cdt, FaceData toReplaceFace) {
  
  // Find the face from FaceData
  CDT::Face_handle f1 = get_face_from_face_data(cdt, toReplaceFace);

  // Initialize the return value
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
  std::list<FaceData> faces;
  faces.push_back(toReplaceFace);
  if (mergable_neigh1) {
    mark_points_to_remove(cdt, e1, neigh1, edges_to_remove, to_remove_points, polygon_cdt, faces);
  }
  if (mergable_neigh2) {
    mark_points_to_remove(cdt, e2, neigh2, edges_to_remove, to_remove_points, polygon_cdt, faces);
  }
  if (mergable_neigh3) {
    mark_points_to_remove(cdt, e3, neigh3, edges_to_remove, to_remove_points, polygon_cdt, faces);
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

  // Save a copy of the cdt
  CDT copy(cdt);

  // Remove the points
  for (const CDT::Vertex_handle v : to_remove_points) {
    cdt.remove_no_flip(v);
  }

  // Add the centroid (or mean point) of the polygon
  cdt.insert_no_flip(centroid);
  cdt.insert_steiner_x_y(centroid.x(), centroid.y());

  // Add all the edges of the merged faces, except from the shared edges as constrains
  std::set<Edge> edges_made_constrained;

  for (FaceData face : faces) {
    for (int i = 0; i < 3; ++i) {

      CDT::Face_handle temp_face = get_face_from_face_data(copy, face);
      Point temp_p1 = temp_face->vertex((i + 1) % 3)->point();
      Point temp_p2 = temp_face->vertex((i + 2) % 3)->point();
      Edge e = std::make_pair(temp_face, i);
      
      if (cdt.is_constrained(e)) { // if the edge is already constrained skip
        continue;
      }
      bool to_remove = false;
      for (const auto& edge : edges_to_remove) {
        if (equal_edges(temp_p1, temp_p2, edge.first, edge.second)) {
          to_remove = true;
          break;
        }
      }
      if (!to_remove) {
        cdt.insert_constraint(temp_p1, temp_p2);
        edges_made_constrained.insert(e);
        break;
      }
    }
  }

  // Remove the shared edges as constraints
  for (const auto& e : edges_made_constrained) {
    cdt.remove_constrained_edge(e.first, e.second);
  }

  ret.obt_count = count_obtuse_triangles(cdt);
  std::list<CDT::Face_handle> affected_faces;
  affected_faces.push_back(f1);
  if (mergable_neigh1) affected_faces.push_back(neigh1);
  if (mergable_neigh2) affected_faces.push_back(neigh2);
  if (mergable_neigh3) affected_faces.push_back(neigh3);
  ret.affected_faces = affected_faces;
  return ret;
}