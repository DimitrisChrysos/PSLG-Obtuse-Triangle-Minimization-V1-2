#include "includes/custom_cdt_class/custom_cdt_class.hpp"
#include "includes/utils/utils.hpp"
#include "includes/read_write_file/read_write_file.hpp"
#include "includes/steiner_methods/steiner_methods.hpp"
#include <random>
#include <chrono>

typedef CGAL::Constrained_triangulation_plus_2<custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag>> CDT;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment;
typedef CDT::Point Point;
typedef CDT::Edge Edge;

using namespace utils;
using namespace steiner_methods;




// Accept or decline something with the given probability
bool accept_or_decline(double prob) {
  std::random_device rd; // Seed for random number generator
  std::mt19937 gen(rd()); // Mersenne Twister random number generator
  std::bernoulli_distribution dist(prob); // Bernoulli distribution
  
  // Generate the random value
  return dist(gen);
}

InsertionMethod choose_steiner_method() {
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

using SteinerMethodObtPoint = obt_point (*)(CDT&, CDT::Face_handle);
using SteinerMethodObtFace = obt_face (*)(CDT&, CDT::Face_handle);

void sim_annealing(CDT& cdt, double a, double b, int L) {
  int steiner_counter = 0;
  double cur_en = a * count_obtuse_triangles(cdt) + b * steiner_counter;
  double T = 1;
  double new_en, de;
  double e = std::exp(1);
  while (T >= 0) {
    
    bool flag = true;

    while (flag) {
      flag = false;

      if (T < 0) {
        break;
      }

      // Iterate the faces of the cdt
      for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); fit++) {

        CDT::Face_handle f = fit;

        if (!is_triangle_inside_region_boundary(f))
          continue;

        if (has_obtuse_angle(f)) {
          
          std::cout<<"1.\n";

          InsertionMethod steiner_method = choose_steiner_method();
          std::cout << "steiner_method: " << static_cast<int>(steiner_method) << std::endl;
          if (steiner_method == InsertionMethod::PROJECTION || 
              steiner_method == InsertionMethod::MIDPOINT || 
              steiner_method == InsertionMethod::CENTROID) {

            SteinerMethodObtPoint method;
            if (steiner_method == InsertionMethod::PROJECTION) method = insert_projection;
            else if (steiner_method == InsertionMethod::MIDPOINT) method = insert_mid;
            else if (steiner_method == InsertionMethod::CENTROID) method = insert_centroid;

            std::cout<<"2.\n";

            CDT copy(cdt);
            obt_point calc_insert_proj = method(copy, f);
            steiner_counter++;
            new_en = a * calc_insert_proj.obt_count + b * steiner_counter;
            de = new_en - cur_en;
            if (de < 0) {
              method(cdt, f);
              cur_en = new_en;
              std::cout<<"3.1.\n";
              flag = true;
              break;
            }
            else {
              double exponent = (-1*de) / T;
              double prob = std::pow(e, exponent);
              bool acc = accept_or_decline(prob);
              std::cout<<"3.2.\n";
              if (acc) {
                method(cdt, f);
                cur_en = new_en;
                flag = true;
                break;
              }
              else {
                steiner_counter--;
              }
            }
          }
          else if (steiner_method == InsertionMethod::CIRCUMCENTER ||
                    steiner_method == InsertionMethod::MERGE_OBTUSE) {

            SteinerMethodObtFace method;
            if (steiner_method == InsertionMethod::CIRCUMCENTER) method = insert_circumcenter;
            else if (steiner_method == InsertionMethod::MERGE_OBTUSE) method = merge_obtuse;

            std::cout<<"4.\n";
            //TODO EDO EIMASTE!!!!
            CDT copy(cdt);
            obt_face temp = method(copy, f);
            std::cout<<"4.123\n";
            if (temp.obt_count == 9999) continue;
            steiner_counter++;
            new_en = a * temp.obt_count + b * steiner_counter;
            de = new_en - cur_en;
            std::cout<<"5.\n";
            if (de < 0) {
              method(cdt, f);
              cur_en = new_en;
              std::cout<<"5.1\n";
              flag = true;
              break;
            }
            else {
              double exponent = (-1*de) / T;
              double prob = std::pow(e, exponent);
              bool acc = accept_or_decline(prob);
              std::cout<<"5.2\n";
              if (acc) {
                std::cout<<"5.2.1\n";
                method(cdt, f);
                cur_en = new_en;
                std::cout<<"5.2.2\n";
                flag = true;
                break;
              }
              else {
                steiner_counter--;
              }
            }
          }
        }
      }

      if (count_obtuse_triangles(cdt) == 0) {
        std::cout << "current obtuse triangles: 0\n" << std::endl;
        return;
      }

      T = T - (double)((double)1 / (double)L);
      std::cout << "After try | obt_triangles: " << count_obtuse_triangles(cdt) << " | steiner_counter: " << steiner_counter << "| T: " << T << std::endl;
    }
    
    // }
  }
}


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

    auto start = std::chrono::high_resolution_clock::now();
    
    // Iterate the faces of the cdt
    for (CDT::Finite_faces_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++) {
      
      // CDT::Face_handle f = fit;
      

      if (!is_triangle_inside_region_boundary(f))
        continue;

      if (has_obtuse_angle(f)) {
        

        CDT copy4(cdt);
        // CDT::Face_handle copy_face1 = find_matching_face(copy4, f);
        obt_face temp_merge_face = merge_obtuse(copy4, f);
        if (best_steiner.obt_count > temp_merge_face.obt_count) {
          merge_face = temp_merge_face;
          best_method = InsertionMethod::MERGE_OBTUSE;
        }

        // CDT::Face_handle copy_face2 = find_matching_face(copy3, f);
        CDT copy3(cdt);
        obt_face temp_circ_face = insert_circumcenter(copy3, f);
        if (best_steiner.obt_count > temp_circ_face.obt_count) {
          circumcenter_face = temp_circ_face;
          best_method = InsertionMethod::CIRCUMCENTER;
        }

        // std::cout << "1. XAAX!\n";

        CDT copy2(cdt);
        obt_point calc_insert_centr = insert_centroid(copy2, f);
        if (best_steiner.obt_count > calc_insert_centr.obt_count) {
          best_steiner = calc_insert_centr;
          best_method = InsertionMethod::CENTROID;
        }

        // std::cout << "2. XAAX!\n";


        CDT copy1(cdt);
        obt_point calc_insert_mid = insert_mid(copy1, f);
        if (best_steiner.obt_count > calc_insert_mid.obt_count) {
          best_steiner = calc_insert_mid;
          best_method = InsertionMethod::MIDPOINT;
        }

        // std::cout << "3. XAAX!\n";


        CDT copy(cdt);
        obt_point calc_insert_proj = insert_projection(copy, f);
        if (best_steiner.obt_count > calc_insert_proj.obt_count) {
          best_steiner = calc_insert_proj;
          best_method = InsertionMethod::PROJECTION;
        }

        // std::cout << "4. XAAX!\n";

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
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    if (elapsed.count() > 1) {
      std::cout << "Time limit reached\n";
      break;
    }
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
    if (!delaunay) {
      local_search(cdt, 30);
      if (count_obtuse_triangles(cdt) == 0) {
        return;
      }
      local_search(cdt, L);
    }
    else {
      local_search(cdt, L);
    }
  }
  else if (method == "sa") {
    auto it = parameters.begin();
    double alpha = it->second;
    it++;
    double beta = it->second;
    it++;
    double L = it->second;
    if (!delaunay) {
      local_search(cdt, 30);
      if (count_obtuse_triangles(cdt) == 0) {
        return;
      }
      sim_annealing(cdt, alpha, beta, L);
    }
    else {
      sim_annealing(cdt, alpha, beta, L);
    }
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
    if (!delaunay) {
      local_search(cdt, 30);
      if (count_obtuse_triangles(cdt) == 0) {
        return;
      }
      // ant method
    }
    else {
      // ant method
    }
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
  boost::property_tree::ptree parameters_for_output = root.get_child("parameters");
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
  // CGAL::draw(cdt);
  handle_methods(cdt, method, parameters, delaunay);

  write_output(cdt, points, method, parameters_for_output, argv[4]);
  // CGAL::draw(cdt);
  
  return 0;
}