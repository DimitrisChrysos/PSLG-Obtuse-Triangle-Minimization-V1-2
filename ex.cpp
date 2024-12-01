#include "includes/custom_cdt_class/custom_cdt_class.hpp"
#include "includes/utils/utils.hpp"
#include "includes/read_write_file/read_write_file.hpp"
#include "includes/steiner_methods/steiner_methods.hpp"
#include <random>

typedef CGAL::Constrained_triangulation_plus_2<custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag>> CDT;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment;
typedef CDT::Point Point;
typedef CDT::Edge Edge;

using namespace utils;
using namespace steiner_methods;



bool accept_or_decline(double prob) {
  std::random_device rd; // Seed for random number generator
  std::mt19937 gen(rd()); // Mersenne Twister random number generator
  std::bernoulli_distribution dist(prob); // Bernoulli distribution
  
  // Generate the random value
  return dist(gen);
}

// int sim_annealing(CDT& cdt, double a, double b, int L) {
//   int steiner_counter;
//   double cur_en = a * count_obtuse_triangles(cdt) + b * steiner_counter;
//   double T = 1;
//   double new_en, de;
//   double e = std::exp(1);
//   while (T >= 0) {
//     bool flag = true;

//     while (flag) {
//       flag = false;

//       // Iterate the faces of the cdt
//       for (CDT::Finite_faces_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++) {

//         if (!is_triangle_inside_region_boundary(f))
//           continue;

//         if (has_obtuse_angle(f)) {
          
//           CDT copy(cdt);
//           obt_point calc_insert_proj = insert_projection(copy, f);
//           steiner_counter++;
//           new_en = a * calc_insert_proj.obt_count + b * steiner_counter;
//           de = new_en - cur_en;
//           if (de < 0) {
//             insert_projection(cdt, f);
//             cur_en = new_en;
//             flag = true;
//             break;
//           }
//           else {
//             double exponent = (-1*de) / T;
//             double prob = std::pow(e, exponent);
//             bool acc = accept_or_decline(prob);
//             if (acc) {
//               insert_projection(cdt, f);
//               cur_en = new_en;
//               flag = true;
//               break;
//             }
//             else {
//               steiner_counter--;
//             }
//           }

//           CDT copy1(cdt);
//           obt_point calc_insert_mid = insert_mid(copy1, f);
//           steiner_counter++;
//           new_en = a * calc_insert_mid.obt_count + b * steiner_counter;
//           de = new_en - cur_en;
//           if (de < 0) {
//             insert_mid(cdt, f);
//             cur_en = new_en;
//             flag = true;
//             break;
//           }
//           else {
//             double exponent = (-1*de) / T;
//             double prob = std::pow(e, exponent);
//             bool acc = accept_or_decline(prob);
//             if (acc) {
//               insert_mid(cdt, f);
//               cur_en = new_en;
//               flag = true;
//               break;
//             }
//             else {
//               steiner_counter--;
//             }
//           }

//           CDT copy2(cdt);
//           obt_point calc_insert_centr = insert_centroid(copy2, f);
//           steiner_counter++;
//           new_en = a * calc_insert_centr.obt_count + b * steiner_counter;
//           de = new_en - cur_en;
//           if (de < 0) {
//             insert_centroid(cdt, f);
//             cur_en = new_en;
//             flag = true;
//             break;
//           }
//           else {
//             double exponent = (-1*de) / T;
//             double prob = std::pow(e, exponent);
//             bool acc = accept_or_decline(prob);
//             if (acc) {
//               insert_centroid(cdt, f);
//               cur_en = new_en;
//               flag = true;
//               break;
//             }
//             else {
//               steiner_counter--;
//             }
//           }

//           CDT copy3(cdt);
//           obt_point calc_insert_circ = insert_circumcenter(copy3, f);
//           if (calc_insert_circ.obt_count == 9999) continue;
//           steiner_counter++;
//           new_en = a * calc_insert_circ.obt_count + b * steiner_counter;
//           de = new_en - cur_en;
//           if (de < 0) {
//             insert_circumcenter(cdt, f);
//             cur_en = new_en;
//             flag = true;
//             break;
//           }
//           else {
//             double exponent = (-1*de) / T;
//             double prob = std::pow(e, exponent);
//             bool acc = accept_or_decline(prob);
//             if (acc) {
//               insert_circumcenter(cdt, f);
//               cur_en = new_en;
//               flag = true;
//               break;
//             }
//             else {
//               steiner_counter--;
//             }
//           }
          

//           CDT copy4(cdt);
//           obt_face temp = merge_obtuse(copy4, f);
//           if (temp.obt_count == 9999) continue;
//           steiner_counter++;
//           new_en = a * temp.obt_count + b * steiner_counter;
//           de = new_en - cur_en;
//           if (de < 0) {
//             merge_obtuse(cdt, f);
//             cur_en = new_en;
//             flag = true;
//             break;
//           }
//           else {
//             double exponent = (-1*de) / T;
//             double prob = std::pow(e, exponent);
//             bool acc = accept_or_decline(prob);
//             if (acc) {
//               merge_obtuse(cdt, f);
//               cur_en = new_en;
//               flag = true;
//               break;
//             }
//             else {
//               steiner_counter--;
//             }
//           }
//         }
//       }
//       T = T - (1 / L);
//     }
//   }
// }


// Local Search for steiner insertions
void local_search(CDT& cdt, int L) {
  int i;
  for (i = 0 ; i < L ; i++) {

    if (count_obtuse_triangles(cdt) == 0) {
      break;
    }

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
    }
    else if (best_method == InsertionMethod::CIRCUMCENTER) {
      insert_circumcenter(cdt, circumcenter_face.face);
    }
    else if (best_method == InsertionMethod::MERGE_OBTUSE) {
      merge_obtuse(cdt, merge_face.face);
    }
    std::cout << "After try " << i << " to insert Steiner | obt_triangles: " << count_obtuse_triangles(cdt) << std::endl;
  }
  std::cout << "After " << i << " Steiner Insertions | obt_triangles: " << count_obtuse_triangles(cdt) << std::endl;
}


// Create the region_boundary_polygon polygon
Polygon_2 make_region_boundary_polygon(std::list<int> region_boundary, std::vector<Point> points) {
  Polygon_2 region_boundary_polygon;
  for (int temp : region_boundary) {
    region_boundary_polygon.push_back(points[temp]);
  }
  return region_boundary_polygon;
}

void handle_methods(CDT& cdt, 
                    std::string method, 
                    std::list<std::pair<std::string, double>> parameters,
                    bool delaunay) {
  if (method == "local") {
    auto it = parameters.begin();
    double L = it->second;
    local_search(cdt, L);
  }
  else if (method == "sa") {
    auto it = parameters.begin();
    double alpha = it->second;
    it++;
    double beta = it->second;
    it++;
    double L = it->second;
    // sim_annealing(cdt, alpha, beta, L);
  }
  else if (method == "ant") {
    auto it = parameters.begin();
    double alpha = it->second;
    it++;
    double beta = it->second;
    it++;
    double xi = it->second;
    it++;
    double psi = it->second;
    it++;
    double lambda = it->second;
    it++;
    double kappa = it->second;
    it++;
    double L = it->second;
    //
  }
}

using namespace read_write_file;

int main(int argc, char *argv[]) {
  if (argc != 5) {
    std::cout << "Wrong number of arguments\n";
    return 1;
  }
  
  namespace pt = boost::property_tree; // namespace alias
  pt::ptree root; // create a root node

  // Choose the file to read
  pt::read_json(argv[2], root); // read the json file
  
  // Read the json file
  std::string instance_uid = get_instance_uid(root);
  int num_points = get_num_points(root);
  std::list<int> points_x = get_points_x(root);
  std::list<int> points_y = get_points_y(root);
  std::list<int> region_boundary = get_region_boundary(root);
  std::string num_constraints = get_num_constraints(root);
  std::list<std::pair<int, int>> additional_constraints = get_additional_constraints(root, region_boundary);
  std::string method = get_method(root);
  std::list<std::pair<std::string, double>> parameters = get_parameters(root);
  bool delaunay = get_delaunay(root);

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
  std::cout << "Starting obtuse counter: " << count_obtuse_triangles(cdt) << std::endl;

  // Insert Steiner points
  CGAL::draw(cdt);
  handle_methods(cdt, method, parameters, delaunay);

  write_output(cdt, points, method, parameters, argv[4]);
  CGAL::draw(cdt);
  
  return 0;
}