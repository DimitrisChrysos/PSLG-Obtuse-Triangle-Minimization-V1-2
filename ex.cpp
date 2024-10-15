#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <string>
#include <list>

int main() {
  
  // Read the json file
  namespace pt = boost::property_tree; // namespace alias
  pt::ptree root; // create a root node
  pt::read_json("input.json", root); // read the json file
  
  // // Print the json file
  // pt::write_json(std::cout, root);

  // Read instance_uid
  std::string instance_uid = root.get<std::string>("instance_uid");
  std::cout << "Instance_uid: " << instance_uid << std::endl;

  // Read num_points
  int num_points = root.get<int>("num_points", 0);
  std::cout << "Num_points: " << num_points << std::endl;

  // Read points_x
  std::list<int> points_x;
  for (pt::ptree::value_type &point_x : root.get_child("points_x")) {
    points_x.push_back(point_x.second.get_value<int>());
  }
  // Print points_x
  std::cout << "\nList of points_x:\n";
  for (const auto &point_x : points_x) {
    std::cout << point_x << " ";
  }
  std::cout << std::endl;

  // Read points_y
  std::list<int> points_y;
  for (pt::ptree::value_type &point_y : root.get_child("points_y")) {
    points_y.push_back(point_y.second.get_value<int>());
  }
  // Print points_y
  std::cout << "\nList of points_y:\n";
  for (const auto &point_y : points_y) {
    std::cout << point_y << " ";
  }
  std::cout << std::endl;

  // Read region_boundary
  std::list<int> region_boundary;
  for (pt::ptree::value_type &temp : root.get_child("region_boundary")) {
    region_boundary.push_back(temp.second.get_value<int>());
  }
  // Print region_boundary
  std::cout << "\nList of region_boundary:\n";
  for (const auto &temp : region_boundary) {
    std::cout << temp << " ";
  }
  std::cout << std::endl;

  // Read num_constraints
  std::string num_constraints = root.get<std::string>("num_constraints");
  std::cout << std::endl << "Num_constraints: " << num_constraints << std::endl;

  // Read the additional_constraints
  std::list<std::pair<int, int>> additional_constraints;
  for (pt::ptree::value_type &row : root.get_child("additional_constraints")) {
    auto it = row.second.begin();
    int first = it->second.get_value<int>();
    ++it;
    int second = it->second.get_value<int>();
    additional_constraints.push_back(std::make_pair(first, second));
  }
  // Print the additional_constraints
  std::cout << "\nAdditional_constraints:\n";
  for (const auto &constraint : additional_constraints) {
    std::cout << constraint.first << " " << constraint.second << std::endl;
  }

  return 0;
}