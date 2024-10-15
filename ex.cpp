#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <string>
#include <list>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

typedef CGAL:: Exact_predicates_inexact_constructions_kernel K;
typedef CGAL:: Exact_predicates_tag Itag;
typedef CGAL:: Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Edge Edge;

CDT make_cdt(
              std::list<int> points_x, 
              std::list<int> points_y, 
              std::list<std::pair<int, int>> additional_constraints
              ) {
  // Initialize the Constrained Delaunay Triangulation (CDT) 
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

  // Draw the triangulation using CGAL's draw function
  CGAL::draw(cdt);
  
  return cdt;
}

int main() {
  
  // Read the json file
  namespace pt = boost::property_tree; // namespace alias
  pt::ptree root; // create a root node
  pt::read_json("input.json", root); // read the json file

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

  // Create the Constrained Delaunay Triangulation (CDT)
  CDT cdt = make_cdt(points_x, points_y, additional_constraints);


  return 0;
}