#include "includes/custom_cdt_class/custom_cdt_class.hpp"
#include "includes/utils/utils.hpp"
#include "includes/read_write_file/read_write_file.hpp"
#include "includes/steiner_methods/steiner_methods.hpp"

typedef CGAL::Constrained_triangulation_plus_2<custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag>> CDT;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment;
typedef CDT::Point Point;
typedef CDT::Edge Edge;

using namespace utils;
using namespace steiner_methods;






// Choose the best method to insert a steiner point
int steiner_insertion(CDT& cdt) {
  int init_obtuse_count = count_obtuse_triangles(cdt);
  Point starting_point;
  CDT::Face_handle starting_face;
  obt_point best_steiner(9999, starting_point);
  obt_face merge_face(9999, starting_face);
  obt_face circumcenter_face(9999, starting_face);
  InsertionMethod best_method = InsertionMethod::NONE;

  // Iterate the faces of the cdt
  for (CDT::Finite_faces_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++) {

    if (!is_triangle_inside_region_boundary(f))
      continue;

    if (has_obtuse_angle(f)) {
      
      CDT copy4(cdt);
      obt_face temp_merge_face = merge_obtuse(copy4, f);
      if (best_steiner.obt_count >= temp_merge_face.obt_count) {
        merge_face = temp_merge_face;
        best_method = InsertionMethod::MERGE_OBTUSE;
      }

      CDT copy3(cdt);
      obt_face temp_circ_face = insert_circumcenter(copy3, f);
      if (best_steiner.obt_count >= temp_circ_face.obt_count) {
        circumcenter_face = temp_circ_face;
        best_method = InsertionMethod::CIRCUMCENTER;
      }

      CDT copy2(cdt);
      obt_point calc_insert_centr = insert_centroid(copy2, f);
      if (best_steiner.obt_count >= calc_insert_centr.obt_count) {
        best_steiner = calc_insert_centr;
        best_method = InsertionMethod::CENTROID;
      }

      CDT copy1(cdt);
      obt_point calc_insert_mid = insert_mid(copy1, f);
      if (best_steiner.obt_count >= calc_insert_mid.obt_count) {
        best_steiner = calc_insert_mid;
        best_method = InsertionMethod::MIDPOINT;
      }

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
    insert_circumcenter(cdt, circumcenter_face.face);
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
  std::cout << "Before steiners: " << count_obtuse_triangles(cdt) << std::endl;
  CGAL::draw(cdt);

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