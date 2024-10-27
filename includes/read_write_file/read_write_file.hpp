#include "../custom_cdt_class/custom_cdt_class.hpp"

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

// Create the output file
void write_output(CDT& cdt, std::vector<Point> points) {
  boost::json::object root1;
  root1["content_type"] = "CG_SHOP_2025_Solution";
  root1["instance_uid"] = "some_instance_uid";

  // JSON array to hold steiner_x values
  boost::json::array steiner_x_json;
  for (const auto& x : cdt.steiner_x) {
    double x1 = static_cast<double>(CGAL::to_double(x));
    // If value is an int
    if (x1 == static_cast<double>(static_cast<int>(x1))) {
      steiner_x_json.push_back(static_cast<int>(x1));
    }
    else {
      const auto exact_coord = CGAL::exact(x);
      std::string frac = exact_coord.get_num().get_str() + "/" + exact_coord.get_den().get_str();
      steiner_x_json.push_back(boost::json::value(frac));
    }
  }

  // Add steiner_x JSON array to root object
  root1["steiner_x"] = steiner_x_json;

  // JSON array to hold steiner_y values
  boost::json::array steiner_y_json;
  for (const auto& y : cdt.steiner_y) {
    double y1 = static_cast<double>(CGAL::to_double(y));
    if (y1 == static_cast<double>(static_cast<int>(y1))) {
      steiner_y_json.push_back(static_cast<int>(y1));
    }
    else {
      const auto exact_coord = CGAL::exact(y);
      std::string frac = exact_coord.get_num().get_str() + "/" + exact_coord.get_den().get_str();
      steiner_y_json.push_back(boost::json::value(frac));
    }
  }

  // Add steiner_x JSON array to root object
  root1["steiner_y"] = steiner_y_json;

  boost::json::array edges_array;

  for (CDT::Finite_edges_iterator e = cdt.finite_edges_begin(); e != cdt.finite_edges_end(); e++) {
    CDT::Face_handle f1 = e->first;  // The face containing the edge
    int edge_index = e->second;        // The local index of the edge within the face

    // Get the two vertices of the edge
    CDT::Vertex_handle v1 = f1->vertex((edge_index + 1) % 3);
    CDT::Vertex_handle v2 = f1->vertex((edge_index + 2) % 3);

    // Access the coordinates of the vertices
    Point p1 = v1->point();
    Point p2 = v2->point();

    int ind_1;
    int ind_2;
    bool found_1 = false;
    bool found_2 = false;

    // Traverse the points and find which 2 indexes are contained in the edge
    for (std::size_t i = 0; i < points.size(); ++i) {
      // If you find the first point
      if (p1.x() == points[i].x() && p1.y() == points[i].y()) {
        ind_1 = i;
        found_1 = true;
      }
      if (p2.x() == points[i].x() && p2.y() == points[i].y()) {
        ind_2 = i;
        found_2 = true;
      }
      if (found_1 && found_2) break;
    }

    // If you haven't found both points of the edge, continue with the Steiner points
    if (!found_1 || !found_2) {
      for (std::size_t i = 0; i < cdt.steiner_x.size(); ++i) {
        if (p1.x() == cdt.steiner_x[i] && p1.y() == cdt.steiner_y[i]) {
          // Continue the indexing after the regular points
          ind_1 = points.size() + i;
          found_1 = true;
        }
        if (p2.x() == cdt.steiner_x[i] && p2.y() == cdt.steiner_y[i]) {
          ind_2 = points.size() + i;
          found_2 = true;
        }
        if (found_1 && found_2) break;
      }
    }

    boost::json::array edge_pair;
    edge_pair.push_back(ind_1);
    edge_pair.push_back(ind_2);

    // Add the edge pair to the edges_array
    edges_array.push_back(edge_pair);



  }
  root1["edges"] = edges_array;

  // Serialize the whole JSON object
  std::string json_string = boost::json::serialize(root1);

  // Make the string pretty and ready to print
  std::string pretty_json;
  int indent = 0;
  char prev;
  bool in_string = false;  // Track if we are in a string
  bool in_steiner = false; // Track if we are printing Steiner points or not

  for (char ch : json_string) {
    if (ch == '"') in_string = !in_string;
    if (!in_string) {
      if (ch == ',') {
        if (prev == '"' || prev == ']') {
          pretty_json += ch;
          if (in_steiner) {
            pretty_json += ' ';
          }
          else {
            pretty_json += '\n';
          }
          prev = ch;
          continue;
        }
      }
      else if (ch == ':') {
        pretty_json += ch;
        pretty_json += ' ';
        prev = ch;
        continue;
      }
      else if (ch == '[') {
        in_steiner = true;
        if (prev == '[') {
          in_steiner = false;
          pretty_json += '\n';
        }
        pretty_json += '\t';
      }
      else if (ch == ']') {
        in_steiner = false;
      }
    }
    pretty_json += ch;
    prev = ch;
  }
  
  // Write to file
  std::ofstream file("output.json");
  file << pretty_json;
  file.close();
}