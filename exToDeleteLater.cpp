// // Παράδειγμα με Πολύγωνο και υπολογισμό επιφάνειας
// #include <CGAL/convex_hull_2.h>
// #include <CGAL/Convex_hull_traits_adapter_2.h>
// #include <CGAL/property_map.h>

// #include "includes/utils/utils.hpp"

// typedef K::Point_2                                          Point_2;
// typedef K::Segment_2                                        Segment_2;
// typedef std::vector<Point_2> Points;

// int main()
// {
//   // create a polygon and put some points in it
//   Polygon_2 p, chp;
//   p.push_back(Point_2(8,0));
//   p.push_back(Point_2(8,4));
//   p.push_back(Point_2(4,0));
//   p.push_back(Point_2(0,4));
//   p.push_back(Point_2(0,0));
//   for(const Point_2& p : p.vertices()){
//     std::cout << p << std::endl;
//   }

//   std::cout << "\n\n";
//   // As the range is not light weight, we have to use a reference
//   const Polygon_2::Vertices& range = p.vertices();
//   Points result;
//   for(auto it = range.begin(); it!= range.end(); ++it){
//     std::cout << *it << std::endl;
//   }

//   std::cout << "\n\n";

//   for(const Segment_2& e  : p.edges()){
//     std::cout << e << std::endl;
//   }
//   std::cout << "\n\n";

//   std::cout << p.area() << std::endl;
//   std::cout << "\n\n";

//   CGAL::convex_hull_2(range.begin(), range.end(), std::back_inserter(result));
//   for(auto it = result.begin(); it!= result.end();++it)
//     chp.push_back(*it);
//   std::cout << chp.area() << std::endl;
//   std::cout << "\n\n";

//   std::cout << utils::simple_or_not(p) << std::endl;
//   std::cout << "\n\n";

//   return EXIT_SUCCESS;
// }


///////////////////////////////////////////////


// // Παράδειγμα Τριγωνοποίησης με περιορισμούς
// #include <CGAL/draw_triangulation_2.h>
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Constrained_Delaunay_triangulation_2.h>
// #include <iostream>

// typedef CGAL:: Exact_predicates_inexact_constructions_kernel K;
// typedef CGAL:: Exact_predicates_tag Itag;
// typedef CGAL:: Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag> CDT;
// typedef CDT::Point Point;
// typedef CDT::Edge Edge;
// int
// main()
// {
//   CDT cdt;
//   std::cout << "Inserting a grid of 5x5 constraints " << std::endl;
//   for (int i = 1; i < 6; ++i) 
//     cdt.insert_constraint( Point(0, i), Point(6,1));
//   for (int j = 1; j < 6; ++j) 
//     cdt.insert_constraint( Point (j, 0), Point(j,6)); 
//   assert(cdt.is_valid());
//   int count = 0;
//   for (const Edge& e : cdt.finite_edges())
//     if (cdt.is_constrained(e))
//       ++count;
//   std::cout << "The number of resulting constrained edges is ";
//   std::cout << count << std::endl;
//   CGAL::draw(cdt);
//   return 0;
// }


///////////////////////////////////////////////


// // Τριγωνοποίηση Delaunay Επίπεδου Γράφου Ευθείων Γραμμών (PSLG)
// #include <CGAL/draw_triangulation_2.h>
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Constrained_Delaunay_triangulation_2.h>
// #include <iostream>

// typedef CGAL:: Exact_predicates_inexact_constructions_kernel K;
// typedef CGAL:: Exact_predicates_tag Itag;
// typedef CGAL:: Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag> CDT;
// typedef CDT::Point Point;
// typedef CDT::Edge Edge;

// int main() 
// {
//   // Initialize the Constrained Delaunay Triangulation (CDT) 
//   CDT cdt;

//   // Define the points from the PSLG (x, y coordinates)
//   std::vector<Point> points = {
//     Point(632, 1588), Point (1330, 1097), Point (3051, 470), Point (5040, 1077), 
//     Point (5883, 2766), Point(8130, 3629), Point (9280, 2836), Point (9613, 4963), 
//     Point(9422, 6363), Point(8996, 7327), Point (8020, 7611), Point(8467, 9720), 
//     Point(6735, 9183), Point(4674, 7865), Point (2519, 7692), Point(973, 9797), 
//     Point (1205, 6005), Point(1929, 5812), Point (3203, 6301), Point (5345, 2923)
//   };

//   // Insert points into the triangulation
//   for (const Point& p : points) {
//     cdt.insert(p);
//   }

//   // Define and add the constrained edges (from additional_constraints)
//   std::vector<std::pair<int, int>> constraints = {
//     {3, 4}, {5, 6}, {9, 10}, {10, 11}, {11, 12}, {12, 13}, {13, 14}, 
//     {14, 15}, {15, 16}, {18, 19}, {19, 0}
//   };

//   // Insert constrained edges based on the provided indices
//   for (const auto& constraint : constraints) {
//     cdt.insert_constraint(points[constraint.first], points[constraint.second]);
//   }

//   // Draw the triangulation using CGAL's draw function
//   CGAL::draw(cdt);
  
//   return 0;
// }


/////////////////////////////////////////////////


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <string>

int main() {
  
  // Read the json file
  namespace pt = boost::property_tree; // namespace alias
  pt::ptree root; // create a root node
  pt::read_json("inputToDeleteLater.json", root); // read the json file

  // Read height
  int height = root.get<int>("height", 0);
  std::cout << "Height: " << height << std::endl;

  // Read message from complex path
  std::string msg = root.get<std::string>("some.complex.path");
  std::cout << "Message: " << msg << std::endl;

  // Read animals
  std::vector< std::pair<std::string, std::string> > animals;
  for (pt::ptree::value_type &animal : root.get_child("animals")) {

    std::string name = animal.first; // get the key
    std::string color = animal.second.data(); // get the value
    animals.push_back(std::make_pair(name, color)); // add to the vector "animals"
  }
  // Print the animals
  std::cout << "\nList of Animals:\n";
  for (const auto &animal : animals) {
      std::cout << "Name: " << animal.first << ", Color: " << animal.second << std::endl;
  }

  // Read fruits
  std::vector<std::string> fruits;
  for (pt::ptree::value_type &fruit : root.get_child("fruits")) {
    fruits.push_back(fruit.second.data());
  }
  // Print the fruits
  std::cout << "\nList of Fruits:\n";
  for (const auto &fruit : fruits) {
    std::cout << fruit << std::endl;
  }

  // Read the matrix
  int matrix[3][3];
  int x = 0;
  for (pt::ptree::value_type &row : root.get_child("matrix")) {
    int y = 0;
    for (pt::ptree::value_type &cell : row.second) {
      matrix[x][y] = cell.second.get_value<int>();
      y++;
    }
    x++;
  }
  // Print the matrix
  std::cout << "\nMatrix:\n";
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      std::cout << matrix[i][j] << " ";
    }
    std::cout << std::endl;
  }

  // Print the json file
  pt::write_json(std::cout, root);




  // Write the json file
  pt::ptree root1;

  // Add new values to the json file
  root1.put("height", height);
  root1.put("some.complex.path", "bonjour");

  // Add animals
  pt::ptree animals_node; // create a new node
  for (const auto &animal : animals) { // add animals as childs
    animals_node.put(animal.first, animal.second);
  }
  root1.add_child("animals", animals_node); // add the node to the root1

  // Add fruits
  pt::ptree fruits_node;
  for (auto &fruit : fruits) {

    // Create a new node for each fruit
    pt::ptree fruit_node;
    fruit_node.put("", fruit);

    // Add the fruit node to the fruits node
    fruits_node.push_back(std::make_pair("", fruit_node));
  }
  root1.add_child("fruits", fruits_node);

  // Add a matrix
  pt::ptree matrix_node;
  for (int i = 0; i < 3; i++) {

    // Create a new node for each row
    pt::ptree row_node;
    for (int j = 0; j < 3; j++) {

      // Create a new node for each cell
      pt::ptree cell_node;

      // Add the value to the cell node
      cell_node.put("", matrix[i][j]);

      // Add the cell node to the row node
      row_node.push_back(std::make_pair("", cell_node));
    }

    // Add the row node to the matrix node
    matrix_node.push_back(std::make_pair("", row_node));
  }

  // Add the matrix node to the root1
  root1.add_child("matrix", matrix_node);

  // Print the json file
  pt::write_json(std::cout, root1);

  // Write the json file
  pt::write_json("output.json", root1);
}