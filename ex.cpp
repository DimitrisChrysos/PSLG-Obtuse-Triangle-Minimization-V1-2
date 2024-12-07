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


void store_face_data(CDT::Face_handle face, std::list<FaceData>& affected_faces) {
  Point p1 = face->vertex(0)->point();
  Point p2 = face->vertex(1)->point();
  Point p3 = face->vertex(2)->point();
  affected_faces.emplace_back(p1, p2, p3);
}

class ant_parameters {
    public:
      double alpha;
      double beta;
      double xi;
      double psi;
      double lambda;
      double kappa;
      double L;

      ant_parameters(double alpha, double beta, double xi, double psi, double lambda, double kappa, double L) {
        this->alpha = alpha;
        this->beta = beta;
        this->xi = xi;
        this->psi = psi;
        this->lambda = lambda;
        this->kappa = kappa;
        this->L = L;
      }
};

class t_sp {
    public:
      double projection;
      double midpoint;
      double centroid;
      double circumcenter;
      double merge_obtuse;

      t_sp(double projection, double midpoint, double centroid, double circumcenter, double merge_obtuse) {
        this->projection = projection;
        this->midpoint = midpoint;
        this->centroid = centroid;
        this->circumcenter = circumcenter;
        this->merge_obtuse = merge_obtuse;
      }
};

class dt {
    public:
      double projection;
      double midpoint;
      double centroid;
      double circumcenter;
      double merge_obtuse;

      dt(double projection, double midpoint, double centroid, double circumcenter, double merge_obtuse) {
        this->projection = projection;
        this->midpoint = midpoint;
        this->centroid = centroid;
        this->circumcenter = circumcenter;
        this->merge_obtuse = merge_obtuse;
      }
};

int ant_id = 0;
class effective_ant {
    public:
      int id;
      int obt_count;
      int steiner_count;
      InsertionMethod sp_method;
      std::list<FaceData> affected_faces;
      Point insrt_pt;
      FaceData face_for_sp_method = FaceData(Point(0,0), Point(0,0), Point(0,0));

      effective_ant(int obt_count, int steiner_count, InsertionMethod sp_method, std::list<FaceData> affected_faces)
        : id(ant_id++), obt_count(obt_count), steiner_count(steiner_count), sp_method(sp_method), affected_faces(affected_faces) {}

      bool operator==(const effective_ant& other) const {
        return id == other.id;
      }
};

void test_affected_faces(std::list<FaceData> affected_faces) {
  for (auto face : affected_faces) {
    std::cout << "p1: " << face.p1 << " p2: " << face.p2 << " p3: " << face.p3 << std::endl;
  }
}

double eucledean_distance(Point p1, Point p2) {
  double p1x = CGAL::to_double(p1.x());
  double p1y = CGAL::to_double(p1.y());
  double p2x = CGAL::to_double(p2.x());
  double p2y = CGAL::to_double(p2.y());
  return std::sqrt(std::pow(p1x - p2x, 2) + std::pow(p1y - p2y, 2));
}

double largest_edge_length(CDT::Face_handle face) {
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

double triangle_height_from_longest_side(CDT::Face_handle face) {
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
double calculate_r_to_h(CDT::Face_handle face) {
  return (double)(( (double)largest_edge_length(face) / (double)2) / (double)triangle_height_from_longest_side(face));
}

bool more_or_equal_to_2_adjacent_obtuse_faces(CDT& cdt, CDT::Face_handle face) {
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



double calculate_posibility(double t, double h, double xi, double psi, double sum) {
  return (double)( std::pow(t, xi) * std::pow(h, psi) ) / sum;
}

InsertionMethod choose_steiner_method(CDT& cdt, CDT::Face_handle face, double k, double xi, double psi, t_sp tsp) {

  // Calculate the radius-to-height ratio
  double r_to_h = calculate_r_to_h(face);

  // Heuristics
  double h_projection = std::max(0.0 , (double)(r_to_h - (double)1) / r_to_h );
  double h_circumcenter = (double) r_to_h / (double)((double)2 + r_to_h);
  double h_midpoint = std::max(0.0, ((double)3 - (double)((double)2 * r_to_h) ) / (double)3 );
  double h_merge_obtuse = more_or_equal_to_2_adjacent_obtuse_faces(cdt, face) ? 1 : 0;

  // Posibilities
  double sum = pow(tsp.projection, xi)*pow(h_projection, psi) + 
                pow(tsp.circumcenter, xi)*pow(h_circumcenter, psi) + 
                pow(tsp.midpoint, xi)*pow(h_midpoint, psi) + 
                pow(tsp.merge_obtuse, xi)*pow(h_merge_obtuse, psi);
  double p_projection = calculate_posibility(tsp.projection, h_projection, xi, psi, sum);
  double p_circumcenter = calculate_posibility(tsp.circumcenter, h_circumcenter, xi, psi, sum);
  double p_midpoint = calculate_posibility(tsp.midpoint, h_midpoint, xi, psi, sum);
  double p_merge_obtuse = calculate_posibility(tsp.merge_obtuse, h_merge_obtuse, xi, psi, sum);

  // Choose the method
  double sum_of_probabilities = p_projection + p_circumcenter + p_midpoint + p_merge_obtuse;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, sum_of_probabilities);
  double random_number = dis(gen);
  // std::cout << "P_proj: " << p_projection << " P_circ: " << p_circumcenter << " P_mid: " << p_midpoint << " P_merge: " << p_merge_obtuse << std::endl;
  if (random_number <= p_projection) {
    // std::cout << "Einai projection" << std::endl;
    return InsertionMethod::PROJECTION;
  }
  else if (random_number <= p_projection + p_circumcenter) {
    // std::cout << "Einai circumcenter" << std::endl;
    return InsertionMethod::CIRCUMCENTER;
  }
  else if (random_number <= p_projection + p_circumcenter + p_midpoint) {
    // std::cout << "Einai midpoint" << std::endl;
    return InsertionMethod::MIDPOINT;
  }
  else {
    // std::cout << "Einai merge obtuse" << std::endl;
    return InsertionMethod::MERGE_OBTUSE;
  }
}

effective_ant improve_trianglulation(CDT& cdt, double k, ant_parameters ant_params, t_sp tsp) {

  // Get a random face
  std::list<CDT::Face_handle> obtuse_faces;
  for (CDT::Finite_faces_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++) {
    if (is_triangle_inside_region_boundary(f) && has_obtuse_angle(f)) {
      obtuse_faces.push_back(f);
    }
  }
  if (obtuse_faces.empty()) { // Fails because the triangluation is already optimal
    return effective_ant(0, 0, InsertionMethod::NONE, std::list<FaceData>());
  }
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, obtuse_faces.size() - 1);
  int random_index = dis(gen);
  auto it = std::next(obtuse_faces.begin(), random_index);
  CDT::Face_handle random_face = *it;

  // Choose a steiner method:
  InsertionMethod steiner_method = choose_steiner_method(cdt, random_face, k, ant_params.xi, ant_params.psi, tsp);
  SteinerMethodObtPoint method_point;
  SteinerMethodObtFace method_face;
  if (steiner_method == InsertionMethod::PROJECTION) method_point = insert_projection;
  else if (steiner_method == InsertionMethod::MIDPOINT) method_point = insert_mid;
  else if (steiner_method == InsertionMethod::CENTROID) method_point = insert_centroid;
  else if (steiner_method == InsertionMethod::CIRCUMCENTER) method_face = insert_circumcenter;
  else if (steiner_method == InsertionMethod::MERGE_OBTUSE) method_face = merge_obtuse;

  // If the choosen method fails, pick the centroid method
  if (steiner_method == InsertionMethod::CIRCUMCENTER || steiner_method == InsertionMethod::MERGE_OBTUSE) {
    CDT copy(cdt);
    FaceData face_data(random_face->vertex(0)->point(), random_face->vertex(1)->point(), random_face->vertex(2)->point());
    obt_face temp = method_face(copy, face_data);
    if (temp.obt_count == 9999) {
      steiner_method = InsertionMethod::CENTROID;
      method_point = insert_centroid;
    }
  }

  // Use the method
  if (steiner_method == InsertionMethod::PROJECTION || 
      steiner_method == InsertionMethod::MIDPOINT || 
      steiner_method == InsertionMethod::CENTROID) {
    CDT save_cdt(cdt);
    obt_point returned = method_point(cdt, random_face);
    CDT::Face_handle located_face = save_cdt.locate(returned.insrt_pt);
    std::list<FaceData> affected_faces;
    store_face_data(located_face, affected_faces);
    effective_ant returned_ant(count_obtuse_triangles(cdt), 0, steiner_method, affected_faces);
    returned_ant.insrt_pt = returned.insrt_pt;
    return returned_ant;
  }
  else if (steiner_method == InsertionMethod::CIRCUMCENTER || 
           steiner_method == InsertionMethod::MERGE_OBTUSE) {
    FaceData face_data_final(random_face->vertex(0)->point(), random_face->vertex(1)->point(), random_face->vertex(2)->point());
    obt_face returned = method_face(cdt, face_data_final);
    std::list<FaceData> affected_faces;
    for (CDT::Face_handle face : returned.affected_faces) {
      store_face_data(face, affected_faces);
    }
    effective_ant returned_ant(count_obtuse_triangles(cdt), 0, steiner_method, affected_faces);
    FaceData face_data(returned.face->vertex(0)->point(), returned.face->vertex(1)->point(), returned.face->vertex(2)->point());
    returned_ant.face_for_sp_method = face_data;
  return returned_ant;
  }

  // Default return, WARNING! Should not be used!
  return effective_ant(9999, 0, InsertionMethod::NONE, std::list<FaceData>());
}


double update_dt(int obt_count, double dt, int steiner_counter, ant_parameters ant_params) {
  dt = dt + ((double)1 / ( (double)1 + ant_params.alpha*(double)obt_count + ant_params.beta*(double)steiner_counter ));
  return dt;
}

bool evaluate_trianguation(CDT& cdt, int before_cycle_obt_count, InsertionMethod steiner_method, dt& Dt, int steiner_counter, ant_parameters ant_params) {
  int obt_count = count_obtuse_triangles(cdt);
  if (obt_count < before_cycle_obt_count) {
    if (steiner_method == InsertionMethod::PROJECTION) Dt.projection = update_dt(obt_count, Dt.projection, steiner_counter, ant_params);
    else if (steiner_method == InsertionMethod::MIDPOINT) Dt.midpoint = update_dt(obt_count, Dt.midpoint, steiner_counter, ant_params);
    else if (steiner_method == InsertionMethod::CENTROID) Dt.centroid = update_dt(obt_count, Dt.centroid, steiner_counter, ant_params);
    else if (steiner_method == InsertionMethod::CIRCUMCENTER) Dt.circumcenter = update_dt(obt_count, Dt.circumcenter, steiner_counter, ant_params);
    else if (steiner_method == InsertionMethod::MERGE_OBTUSE) Dt.merge_obtuse = update_dt(obt_count, Dt.merge_obtuse, steiner_counter, ant_params);
    return true;
  }
  return false;
}

void update_pheromones(t_sp& tsp, ant_parameters ant_params, dt Dt) {
  tsp.projection = ((double)1 - ant_params.lambda)*tsp.projection + Dt.projection;
  tsp.midpoint = ((double)1 - ant_params.lambda)*tsp.midpoint + Dt.midpoint;
  tsp.centroid = ((double)1 - ant_params.lambda)*tsp.centroid + Dt.centroid;
  tsp.circumcenter = ((double)1 - ant_params.lambda)*tsp.circumcenter + Dt.circumcenter;
  tsp.merge_obtuse = ((double)1 - ant_params.lambda)*tsp.merge_obtuse + Dt.merge_obtuse;
}


// Handles conflicts, returns true if it encounters a conflict
bool handle_conflicts(CDT &cdt, std::list<effective_ant>& effective_ants, effective_ant& new_ant) {
  for (effective_ant& ant : effective_ants) {
    for (FaceData face1 : ant.affected_faces) {
      for (FaceData face2 : new_ant.affected_faces) {
        if (same_faces(face1, face2)) {
          // std::cout << "Vrika Conflict" << std::endl;
          // std::cout << "face1.p1: " << face1.p1 << " face1.p2: " << face1.p2 << " face1.p3: " << face1.p3 << std::endl;
          // std::cout << "face2.p1: " << face2.p1 << " face2.p2: " << face2.p2 << " face2.p3: " << face2.p3 << std::endl;
          if (ant.obt_count > new_ant.obt_count) {
            effective_ants.remove(ant);
            effective_ants.push_back(new_ant);
          }
          return true;
        }
      }
    }
  }
  new_ant.steiner_count += 1;
  effective_ants.push_back(new_ant);
  return false;
}

void use_triangulation_ants(CDT& cdt, std::list<effective_ant>& ants) {
  for (effective_ant& ant : ants) {
    InsertionMethod method = ant.sp_method;
    if (method == InsertionMethod::PROJECTION || 
        method == InsertionMethod::MIDPOINT || 
        method == InsertionMethod::CENTROID) {
      cdt.insert_no_flip(ant.insrt_pt);
      cdt.insert_steiner_x_y(ant.insrt_pt.x(), ant.insrt_pt.y());
    }
    else if (method == InsertionMethod::CIRCUMCENTER) {
      std::cout << "\n1.circum\n";
      insert_circumcenter(cdt, ant.face_for_sp_method); // very likely to fail, because of face not existing
      std::cout << "2.circum\n";
    }
    else if (method == InsertionMethod::MERGE_OBTUSE) {
      std::cout << "1.merge\n";
      merge_obtuse(cdt, ant.face_for_sp_method); // very likely to fail, because of face not existing
      std::cout << "2.merge\n";
    }
  }
}

// Returns true if the triangulation is optimal
// Check to find the best triangulation
bool save_best_triangulation(CDT& cdt, std::list<effective_ant>& effective_ants, 
                              ant_parameters ant_params, double& starting_energy,
                              std::list<effective_ant>& best_triangulation_ants) {
  CDT copy(cdt);
  use_triangulation_ants(copy, effective_ants);
  std::cout << "|| For this cycle obtuse: " << count_obtuse_triangles(copy) << " | Steiner: " << effective_ants.size() << std::endl;
  if (count_obtuse_triangles(copy) == 0) {
    best_triangulation_ants = effective_ants;
    return true;
  }
  int steiner_counter = effective_ants.size();
  double current_energy = ant_params.alpha * count_obtuse_triangles(copy) + ant_params.beta * steiner_counter;
  double de = current_energy - starting_energy;
  if (de < 0) {
    best_triangulation_ants = effective_ants;
    starting_energy = current_energy;
  }
  return false;
}

void ant_colony_optimization(CDT& cdt, ant_parameters ant_params) {

  std::cout<< "L: " << ant_params.L << " | Kappa: " << ant_params.kappa << std::endl;

  if (count_obtuse_triangles(cdt) == 0) return;

  t_sp tsp(0.5, 0.5, 0.5, 0.5, 0.5);
  std::list<effective_ant> best_triangulation_ants;
  double starting_energy = ant_params.alpha * count_obtuse_triangles(cdt) + ant_params.beta * 0;
  for (int c = 1 ; c <= ant_params.L ; c++) {  // c -> cycle

    dt Dt(0, 0, 0, 0, 0);
    int before_cycle_obt_count = count_obtuse_triangles(cdt);
    std::list<effective_ant> effective_ants;

    for (int k = 1 ; k <= ant_params.kappa ; k++) { // k -> ant

      CDT copy(cdt);

      // Improve the triangulation
      effective_ant new_ant = improve_trianglulation(copy, k, ant_params, tsp);
      if (new_ant.obt_count == 0) { // If the triangulation is optimal
        break;
      }
      new_ant.steiner_count = effective_ants.size();

      // Handle conflicts and Evaluate the triangulation
      handle_conflicts(copy, effective_ants, new_ant);
      evaluate_trianguation(copy, before_cycle_obt_count, new_ant.sp_method, Dt, effective_ants.size(), ant_params);
    }

    std::cout << "Cycle: " << c << " finished with -> " << "|| effective ants: " << effective_ants.size();

    // SaveBestTriangulation(c)
    if (save_best_triangulation(cdt, effective_ants, ant_params, starting_energy, best_triangulation_ants)) {
      break; // Break because the triangulation is optimal
    }

    // UpdatePheromones(c)
    update_pheromones(tsp, ant_params, Dt);
  }

  // Use the best triangulation to the starting cdt
  use_triangulation_ants(cdt, best_triangulation_ants);
  std::cout << "\nFinal -> Obtuse Triangles: " << count_obtuse_triangles(cdt) << " || Steiner Points: " << best_triangulation_ants.size() << std::endl;
}

// Accept or decline something with the given probability
bool accept_or_decline(double prob) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::bernoulli_distribution dist(prob);
  return dist(gen); // Generate the random value
}

void sim_annealing(CDT& cdt, double a, double b, int L) {
  int steiner_counter = 0;
  double cur_en = a * count_obtuse_triangles(cdt) + b * steiner_counter;
  double T = 1;
  double new_en, de;
  double e = std::exp(1);
  int tries = 0;
  while (T >= 0) {
    
    tries++;
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
          
          // std::cout<<"1.\n";

          InsertionMethod steiner_method = choose_random_steiner_method();
          // std::cout << "steiner_method: " << static_cast<int>(steiner_method) << std::endl;
          if (steiner_method == InsertionMethod::PROJECTION || 
              steiner_method == InsertionMethod::MIDPOINT || 
              steiner_method == InsertionMethod::CENTROID) {

            SteinerMethodObtPoint method;
            if (steiner_method == InsertionMethod::PROJECTION) method = insert_projection;
            else if (steiner_method == InsertionMethod::MIDPOINT) method = insert_mid;
            else if (steiner_method == InsertionMethod::CENTROID) method = insert_centroid;

            // std::cout<<"2.\n";

            CDT copy(cdt);
            obt_point calc_insert_proj = method(copy, f);
            steiner_counter++;
            new_en = a * calc_insert_proj.obt_count + b * steiner_counter;
            de = new_en - cur_en;
            if (de < 0) {
              method(cdt, f);
              cur_en = new_en;
              // std::cout<<"3.1.\n";
              flag = true;
              break;
            }
            else {
              double exponent = (-1*de) / T;
              double prob = std::pow(e, exponent);
              bool acc = accept_or_decline(prob);
              // std::cout<<"3.2.\n";
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

            // std::cout<<"4.\n";
            //TODO EDO EIMASTE!!!!
            CDT copy(cdt);
            FaceData face_data(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
            obt_face temp = method(copy, face_data);
            // std::cout<<"4.123\n";
            if (temp.obt_count == 9999) continue;
            steiner_counter++;
            new_en = a * temp.obt_count + b * steiner_counter;
            de = new_en - cur_en;
            // std::cout<<"5.\n";
            if (de < 0) {
              FaceData face_data(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
              method(cdt, face_data);
              cur_en = new_en;
              // std::cout<<"5.1\n";
              flag = true;
              break;
            }
            else {
              double exponent = (-1*de) / T;
              double prob = std::pow(e, exponent);
              bool acc = accept_or_decline(prob);
              // std::cout<<"5.2\n";
              if (acc) {
                // std::cout<<"5.2.1\n";
                FaceData face_data(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
                method(cdt, face_data);
                cur_en = new_en;
                // std::cout<<"5.2.2\n";
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
      std::cout << "After sa try " << tries << " -> obt_triangles: " << count_obtuse_triangles(cdt) << " | steiner_counter: " << steiner_counter << "| T: " << T << std::endl;
    }
    
    // }
  }
}


// Local Search for steiner insertions
// Also the method used in the first assignment
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
    FaceData toReplaceFace(Point(0,0), Point(0,0), Point(0,0));

    auto start = std::chrono::high_resolution_clock::now();
    
    // Iterate the faces of the cdt
    for (CDT::Finite_faces_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++) {
      
      // CDT::Face_handle f = fit;
      
      if (!is_triangle_inside_region_boundary(f))
        continue;

      if (has_obtuse_angle(f)) {
        
        CDT copy(cdt);
        obt_point calc_insert_proj = insert_projection(copy, f);
        if (best_steiner.obt_count > calc_insert_proj.obt_count) {
          best_steiner = calc_insert_proj;
          best_method = InsertionMethod::PROJECTION;
        }

        CDT copy1(cdt);
        obt_point calc_insert_mid = insert_mid(copy1, f);
        if (best_steiner.obt_count > calc_insert_mid.obt_count) {
          best_steiner = calc_insert_mid;
          best_method = InsertionMethod::MIDPOINT;
        }

        CDT copy2(cdt);
        obt_point calc_insert_centr = insert_centroid(copy2, f);
        if (best_steiner.obt_count > calc_insert_centr.obt_count) {
          best_steiner = calc_insert_centr;
          best_method = InsertionMethod::CENTROID;
        }

        CDT copy3(cdt);
        FaceData face_data(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
        obt_face temp_circ_face = insert_circumcenter(copy3, face_data);
        if (best_steiner.obt_count > temp_circ_face.obt_count) {
          circumcenter_face = temp_circ_face;
          toReplaceFace.p1 = f->vertex(0)->point();
          toReplaceFace.p2 = f->vertex(1)->point();
          toReplaceFace.p3 = f->vertex(2)->point();
          best_method = InsertionMethod::CIRCUMCENTER;
        }

        CDT copy4(cdt);
        FaceData face_data2(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
        obt_face temp_merge_face = merge_obtuse(copy4, face_data2);
        if (best_steiner.obt_count > temp_merge_face.obt_count) {
          merge_face = temp_merge_face;
          toReplaceFace.p1 = f->vertex(0)->point();
          toReplaceFace.p2 = f->vertex(1)->point();
          toReplaceFace.p3 = f->vertex(2)->point();
          best_method = InsertionMethod::MERGE_OBTUSE;
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
      insert_circumcenter(cdt, toReplaceFace);
    }
    else if (best_method == InsertionMethod::MERGE_OBTUSE) {
      merge_obtuse(cdt, toReplaceFace);
    }
    std::cout << "After local search try " << i << "-> obt_triangles: " << count_obtuse_triangles(cdt) << std::endl;
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
    }
    local_search(cdt, L);
  }
  else if (method == "sa") {
    auto it = parameters.begin();
    double alpha = it->second;
    it++;
    double beta = it->second;
    it++;
    double L = it->second;
    if (!delaunay) {
      local_search(cdt, 70);
      if (count_obtuse_triangles(cdt) == 0) {
        return;
      }
    }
    sim_annealing(cdt, alpha, beta, L);
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
    ant_parameters ant_params(alpha, beta, xi, psi, lambda, kappa, L);
    if (!delaunay) {
      local_search(cdt, 30);
      if (count_obtuse_triangles(cdt) == 0) {
        return;
      }
    }
    ant_colony_optimization(cdt, ant_params);
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
  CGAL::draw(cdt);
  handle_methods(cdt, method, parameters, delaunay);

  write_output(cdt, points, method, parameters_for_output, argv[4]);
  CGAL::draw(cdt);
  
  return 0;
}