#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <string>
#include <list>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS // optional in native ubuntu, removes a warning in wsl
#ifndef CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H

#define CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H



#include <CGAL/Constrained_Delaunay_triangulation_2.h>



template <class Gt, class Tds = CGAL::Default, class Itag = CGAL::Default>

class Custom_Constrained_Delaunay_triangulation_2

    : public CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag> {

public:

    using Base = CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>;

    using typename Base::Face_handle;

    using typename Base::Point;

    using typename Base::Vertex_handle;

    using typename Base::Locate_type;



    // Constructors

    Custom_Constrained_Delaunay_triangulation_2(const Gt& gt = Gt())

        : Base(gt) {}



    Custom_Constrained_Delaunay_triangulation_2(typename Base::List_constraints& lc, const Gt& gt = Gt())

        : Base(lc, gt) {}



    template <class InputIterator>

    Custom_Constrained_Delaunay_triangulation_2(InputIterator it, InputIterator last, const Gt& gt = Gt())

        : Base(it, last, gt) {}



    // New insert method without flips

    Vertex_handle insert_no_flip(const Point& a, Face_handle start = Face_handle()) {

        // Call Ctr::insert without flip_around

        Vertex_handle va = this->Base::Ctr::insert(a, start); // Directly call Ctr::insert from the base

        return va;

    }



    // Another insert method with known location

    Vertex_handle insert_no_flip(const Point& a, Locate_type lt, Face_handle loc, int li) {

        Vertex_handle va = this->Base::Ctr::insert(a, lt, loc, li); // Directly call Ctr::insert from the base

        return va;

    }

    bool my_is_flippable(const typename Base::Edge& e) {
      Face_handle f1 = e.first; 
      int i = e.second; 
      Face_handle f2 = f1->neighbor(i);

      // Check if the edge is a boundary or constrained to mark it as not flippable
      // //
      // // Check if v1 and v2 are valid
      // typename Base::Vertex_handle v1 = f1->vertex((i+1)%3);
      // typename Base::Vertex_handle v2 = f1->vertex((i+2)%3);
      // std::cout << "Edge points: (" << v1->point() << ") - (" << v2->point() << ")" << std::endl;
      // //
      if (this->is_infinite(f1) || this->is_infinite(f2)) {
        return false;
      }
      else if (this->is_constrained(e)) {
        return false;
      }

      return true;
    }
};

#endif // CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H


typedef CGAL:: Exact_predicates_exact_constructions_kernel K;
typedef CGAL:: Exact_predicates_tag Itag;
typedef Custom_Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Edge Edge;

bool has_obtuse_angle(CDT::Face_handle face) {
  // Get the vertices of the triangle
  Point p1 = face->vertex(0)->point();
  Point p2 = face->vertex(1)->point();
  Point p3 = face->vertex(2)->point();

  // Check if any angle of the triangle is obtuse
    if (CGAL::angle(p1,p2,p3) == CGAL::OBTUSE ||
        CGAL::angle(p2,p3,p1) == CGAL::OBTUSE ||
        CGAL::angle(p3,p1,p2) == CGAL::OBTUSE) {
      return true;
    }
  return false;
}

int count_obtuse_triangles(CDT cdt) {
  int count = 0;
  for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); fit++) {
    CDT::Face_handle face = fit;
    
    if (has_obtuse_angle(face)) {
      count++;
    }
  }
  return count;
}

// void make_flips(CDT& cdt) {
//   int count = 0;
//   int save_count;
//   bool flag = true;
//   while (flag) {
//     std::cout << "Irtha edw geiaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n";
//     count = 0;
//     save_count = count;
//     CDT copy(cdt);
//     for (CDT::Finite_edges_iterator e = cdt.finite_edges_begin(); e != cdt.finite_edges_end(); e++) {

//       std::cout << "Now, we have " << count_obtuse_triangles(cdt) << "obtuse triangles\n";
//       // if (count > save_count) break;

//       // Get the 2 triangles
//       CDT::Face_handle f1 = e->first; // The face of the edge
//       int i = e->second; // The index of the edge in the face
//       CDT::Face_handle f2 = f1->neighbor(i); // The face of the neighbor of the edge

//       // Check if the edge is flippable
//       if (!cdt.my_is_flippable(*e)) continue;

//       // Check if the triangles formed by the edge have obtuse angles
//       int obt = 0;
//       if (has_obtuse_angle(f1)) obt++;
//       if (has_obtuse_angle(f2)) obt++;

//       // If triangles have obtuse angles, make the flip
//       if (obt) {
//         // Make the flip
//         cdt.tds().flip(f1, i);
//         count++;

//         // Check if the triangles formed by the edge have obtuse angles after the flip
//         int obt2 = 0;
//         if (has_obtuse_angle(f1)) obt2++;
//         if (has_obtuse_angle(f2)) obt2++;

//         // If the number of obtuse angles is the same or more, undo the flip
//         if (obt2 >= obt) {
//           // std::cout << "obt2: " << obt2 << " | obt: " << obt << std::endl;
//           cdt = copy;
//           count--;
//         }
//         else {
//           std::cout << "obt2: " << obt2 << " | obt: " << obt << std::endl;
//           break;
//           // std::cout << "Now, we have " << count_obtuse_triangles(cdt) << "obtuse triangles\n";
//         }
//       }
//       std::cout << "LOCA LOCA LOCAAAAAA\n";
//     }
//     // std::cout << "Made " << count << " total successful flips" << std::endl;
//     if (count == save_count) flag = false;
//   }
// }

bool test_the_flip(CDT& cdt, Point v1, Point v2) {
  CDT copy(cdt);
  for (CDT::Finite_edges_iterator e = copy.finite_edges_begin(); e != copy.finite_edges_end(); e++) {

    // Get the edge points
    CDT::Face_handle f1 = e->first; 
    int i = e->second; 
    CDT::Face_handle f2 = f1->neighbor(i);
    auto p1 = f1->vertex((i+1)%3);
    auto p2 = f1->vertex((i+2)%3);
    Point p1_point = p1->point();
    Point p2_point = p2->point();

    if (v1.x() == p1_point.x() && v1.y() == p1_point.y() && v2.x() == p2_point.x() && v2.y() == p2_point.y()) {
      std::cout << "Found the edge: (" << p1->point() << ") - (" << p2->point() << ")" << std::endl;
    }

    else 
      continue;

    // Check if the edge is flippable
    if (!copy.my_is_flippable(*e)) return false;

    // Check if the triangles formed by the edge have obtuse angles
    int obt = 0;
    if (has_obtuse_angle(f1)) obt++;
    if (has_obtuse_angle(f2)) obt++;

    // If triangles have obtuse angles, make the flip
    if (obt) {
      // Make the flip
      copy.tds().flip(f1, i);

      // Check if the triangles formed by the edge have obtuse angles after the flip
      int obt2 = 0;
      if (has_obtuse_angle(f1)) obt2++;
      if (has_obtuse_angle(f2)) obt2++;

      // If the number of obtuse angles after the flip are less, return true
      if (obt2 < obt) {
        return true;
      }
    }

    return false;
  }

  return false;
}

void make_flips(CDT& cdt) {
  int count = 0;
  for (CDT::Finite_edges_iterator e = cdt.finite_edges_begin(); e != cdt.finite_edges_end(); e++) {
    
    // Get the edge points
    CDT::Face_handle f1 = e->first; 
    int i = e->second; 
    CDT::Face_handle f2 = f1->neighbor(i);
    auto v1 = f1->vertex((i+1)%3);
    auto v2 = f1->vertex((i+2)%3);
    std::cout << "\nEdge points: (" << v1->point() << ") - (" << v2->point() << ")" << std::endl;

    // Test is the flip possible or if it is worth doing
    bool do_flip = test_the_flip(cdt, v1->point(), v2->point());

    std::cout << "Now, we have " << count_obtuse_triangles(cdt) << " obtuse triangles\n";

    // If the flip is possible and worth it, do it
    if (do_flip) {
      cdt.tds().flip(f1, i);
      count++;
    }
  }
  std::cout << "Made " << count << " total successful flips" << std::endl;
}

// void make_flips(CDT& cdt) {
//   int count = 0;
//   CDT copy = cdt;
//   for (CDT::Finite_edges_iterator e = cdt.finite_edges_begin(); e != cdt.finite_edges_end(); e++) {

//     std::cout << "Now, we have " << count_obtuse_triangles(cdt) << "obtuse triangles\n";

//     // Get the 2 triangles
//     CDT::Face_handle f1 = e->first; // The face of the edge
//     int i = e->second; // The index of the edge in the face
//     CDT::Face_handle f2 = f1->neighbor(i); // The face of the neighbor of the edge

//     // Check if the edge is flippable
//     if (!cdt.my_is_flippable(*e)) continue;

//     // Check if the triangles formed by the edge have obtuse angles
//     int obt = 0;
//     if (has_obtuse_angle(f1)) obt++;
//     if (has_obtuse_angle(f2)) obt++;

//     // If triangles have obtuse angles, make the flip
//     if (obt) {
//       // Make the flip
//       cdt.tds().flip(f1, i);
//       count++;

//       // Check if the triangles formed by the edge have obtuse angles after the flip
//       int obt2 = 0;
//       if (has_obtuse_angle(f1)) obt2++;
//       if (has_obtuse_angle(f2)) obt2++;

//       // If the number of obtuse angles is the same or more, undo the flip
//       if (obt2 >= obt) {
//         // std::cout << "obt2: " << obt2 << " | obt: " << obt << std::endl;
//         cdt = copy;
//         count--;
//       }
//       else {
//         std::cout << "obt2: " << obt2 << " | obt: " << obt << std::endl;
//         // std::cout << "Now, we have " << count_obtuse_triangles(cdt) << "obtuse triangles\n";
//       }
//     }
//   }
//   std::cout << "Made " << count << " total successful flips" << std::endl;
// }

void Steiner_insertion(CDT& cdt) {
  int obt_count = count_obtuse_triangles(cdt);
  std::cout << "Initial obtuse count: " << obt_count << std::endl;
  CDT copy = cdt;
  for (const Edge& e : cdt.finite_edges()) {
    CDT::Face_handle f1 = e.first; // The face of the edge
    if (has_obtuse_angle(f1)) {
      Point a = f1->vertex(0)->point();
      Point b = f1->vertex(1)->point();
      Point c = f1->vertex(2)->point();
      
      std::cout << "before\n";

      std::cout << "a.x -> " << a.x() << "| a.y -> " << a.y() << std::endl;
      std::cout << "b.x -> " << b.x() << "| b.y -> " << b.y() << std::endl;
      std::cout << "c.x -> " << c.x() << "| c.y -> " << c.y() << std::endl;

      Point pericenter = CGAL::circumcenter(a, b, c);
      
      std::cout << "pericenter.x -> " << pericenter.x() << "| pericenter.y -> " << pericenter.y() << std::endl;

      std::cout << "after\n";

      cdt.insert_no_flip(pericenter);
      break;
    }
  }
  int new_obt_count = count_obtuse_triangles(cdt);
  if (new_obt_count >= obt_count) {
    cdt = copy;
  }
}


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

int main() {
  
  // Read the json file
  namespace pt = boost::property_tree; // namespace alias
  pt::ptree root; // create a root node
  pt::read_json("input.json", root); // read the json file
  std::string instance_uid = get_instance_uid(root);
  int num_points = get_num_points(root);
  std::list<int> points_x = get_points_x(root);
  std::list<int> points_y = get_points_y(root);
  std::list<int> region_boundary = get_region_boundary(root);
  std::string num_constraints = get_num_constraints(root);
  std::list<std::pair<int, int>> additional_constraints = get_additional_constraints(root, region_boundary);


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

  

  // // Print all edges
  // for (CDT::Finite_edges_iterator eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); eit++) {
  //   std::cout << "Edge: " << eit->first->vertex((eit->second+1)%3)->point() << " - " << eit->first->vertex((eit->second+2)%3)->point() << std::endl;
  // }

  // Count the obtuse triangles
  int obtuse_angles = count_obtuse_triangles(cdt);
  std::cout << "Number of obtuse triangles before the flips: " << obtuse_angles << std::endl;
  
  // Make flips
  make_flips(cdt);
  // CGAL::draw(cdt);

  // Insert Steiner points
  // CGAL::draw(cdt);
  // Steiner_insertion(cdt);

  // Count the obtuse triangles
  obtuse_angles = count_obtuse_triangles(cdt);
  std::cout << "Number of obtuse triangles after the flips: " << obtuse_angles << std::endl;

  // Draw the triangulation using CGAL's draw function
  CGAL::draw(cdt);

  return 0;
}