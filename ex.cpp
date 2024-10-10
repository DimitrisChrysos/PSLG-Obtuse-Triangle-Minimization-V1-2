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
//   // As the range is not light weight, we have to use a reference
//   const Polygon_2::Vertices& range = p.vertices();
//   Points result;
//   for(auto it = range.begin(); it!= range.end(); ++it){
//     std::cout << *it << std::endl;
//   }
//   for(const Segment_2& e  : p.edges()){
//     std::cout << e << std::endl;
//   }
//   std::cout << p.area() << std::endl;

//   CGAL::convex_hull_2(range.begin(), range.end(), std::back_inserter(result));
//   for(auto it = result.begin(); it!= result.end();++it)
//     chp.push_back(*it);
//   std::cout << chp.area() << std::endl;

//   std::cout << utils::simple_or_not(p) << std::endl;

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


// Τριγωνοποίηση Delaynay Επίπεδου Γράφου Ευθείων Γραμμών (PSLG)
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <iostream>

typedef CGAL:: Exact_predicates_inexact_constructions_kernel K;
typedef CGAL:: Exact_predicates_tag Itag;
typedef CGAL:: Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Edge Edge;

int main() 
{
  // Initialize the Constrained Delaunay Triangulation (CDT) 
  CDT cdt;

  // Define the points from the PSLG (x, y coordinates)
  std::vector<Point> points = {
    Point(632, 1588), Point (1330, 1097), Point (3051, 470), Point (5040, 1077), 
    Point (5883, 2766), Point(8130, 3629), Point (9280, 2836), Point (9613, 4963), 
    Point(9422, 6363), Point(8996, 7327), Point (8020, 7611), Point(8467, 9720), 
    Point(6735, 9183), Point(4674, 7865), Point (2519, 7692), Point(973, 9797), 
    Point (1205, 6005), Point(1929, 5812), Point (3203, 6301), Point (5345, 2923)
  };

  // Insert points into the triangulation
  for (const Point& p : points) {
    cdt.insert(p);
  }

  // Define and add the constrained edges (from additional_constraints)
  std::vector<std::pair<int, int>> constraints = {
    {3, 4}, {5, 6}, {9, 10}, {10, 11}, {11, 12}, {12, 13}, {13, 14}, 
    {14, 15}, {15, 16}, {18, 19}, {19, 0}
  };

  // Insert constrained edges based on the provided indices
  for (const auto& constraint : constraints) {
    cdt.insert_constraint(points[constraint.first], points[constraint.second]);
  }

  // Draw the triangulation using CGAL's draw function
  CGAL::draw(cdt);
  
  return 0;
}