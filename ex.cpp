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

};



#endif // CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H




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
typedef Custom_Constrained_Delaunay_triangulation_2<K, CGAL:: Default, Itag> CDT;
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
  
  return cdt;
}

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

void make_flips(CDT& cdt) {
  int count=0;
  for (const Edge& e : cdt.finite_edges()) {
    CDT::Face_handle f1 = e.first; // The face of the edge
    int i = e.second; // The index of the edge in the face
    CDT::Face_handle f2 = f1->neighbor(i); // The face of the neighbor of the edge

    std::cout << "Coords of f1" << std::endl;

    for (int i = 0; i < 3; ++i) {
        Point p = f1->vertex(i)->point();
        std::cout << "Vertex " << i << ": (" << p.x() << ", " << p.y() << ")\n";
    }

    std::cout << "Coords of f2" << std::endl;

    for (int i = 0; i < 3; ++i) {
        Point p = f2->vertex(i)->point();
        std::cout << "Vertex " << i << ": (" << p.x() << ", " << p.y() << ")\n";
    }

    // Print face handles
    std::cout << "Face handle f1: " << &(*f1) << std::endl;
    std::cout << "Face handle f2: " << &(*f2) << std::endl;

    // // Check if the edge in on the boundary
    // if (cdt.is_infinite(f1) || cdt.is_infinite(f2)) {
    //   continue;
    // }

    // // Check if the edge is constrained
    // if (cdt.is_constrained(e)) {
    //   continue;
    // }


    CDT::Face_handle ni = f1->neighbor(i);
    if (cdt.is_infinite(f1) || cdt.is_infinite(ni)) continue;
    if (f1->is_constrained(i)) continue;

    // Check if the triangles formed by the edge have obtuse angles
    if (has_obtuse_angle(f1) || has_obtuse_angle(f2)) {
      std::cout << count << " Geia sou!\n";
      // if (count == 17) {
      // }
      CGAL::draw(cdt);
      cdt.tds().flip(f1, i);  // Perform a flip
      std::cout << count << " Bye!\n";
    }
    count++;


  }
}


    // if (!cdt.is_flipable(f1, i)) {
    //   continue;
    // }
    // if (has_obtuse_angle(f1) || has_obtuse_angle(f2)) {
    //   count++;
    //   std::cout << "I made " << count << " flips\n";
    //   cdt.tds().flip(f1, i);  // Perform a flip
    // }


int main() {
  
  // Read the json file
  namespace pt = boost::property_tree; // namespace alias
  pt::ptree root; // create a root node
  pt::read_json("input.json", root); // read the json file

  // Read instance_uid
  std::string instance_uid = root.get<std::string>("instance_uid");
  // std::cout << "Instance_uid: " << instance_uid << std::endl;

  // Read num_points
  int num_points = root.get<int>("num_points", 0);
  // std::cout << "Num_points: " << num_points << std::endl;

  // Read points_x
  std::list<int> points_x;
  for (pt::ptree::value_type &point_x : root.get_child("points_x")) {
    points_x.push_back(point_x.second.get_value<int>());
  }
  // // Print points_x
  // std::cout << "\nList of points_x:\n";
  // for (const auto &point_x : points_x) {
  //   std::cout << point_x << " ";
  // }
  // std::cout << std::endl;

  // Read points_y
  std::list<int> points_y;
  for (pt::ptree::value_type &point_y : root.get_child("points_y")) {
    points_y.push_back(point_y.second.get_value<int>());
  }
  // // Print points_y
  // std::cout << "\nList of points_y:\n";
  // for (const auto &point_y : points_y) {
  //   std::cout << point_y << " ";
  // }
  // std::cout << std::endl;

  // Read region_boundary
  std::list<int> region_boundary;
  for (pt::ptree::value_type &temp : root.get_child("region_boundary")) {
    region_boundary.push_back(temp.second.get_value<int>());
  }
  // // Print region_boundary
  // std::cout << "\nList of region_boundary:\n";
  // for (const auto &temp : region_boundary) {
  //   std::cout << temp << " ";
  // }
  // std::cout << std::endl;

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
  // std::cout << "\nAdditional_constraints:\n";
  // for (const auto &constraint : additional_constraints) {
  //   std::cout << constraint.first << " " << constraint.second << std::endl;
  // }

  // Create the Constrained Delaunay Triangulation (CDT)
  CDT cdt = make_cdt(points_x, points_y, additional_constraints);

  // Count the obtuse triangles
  int obtuse_angles = count_obtuse_triangles(cdt);
  std::cout << "\nNumber of obtuse triangles: " << obtuse_angles << std::endl;
  
  // Make flips
  make_flips(cdt);

  // Count the obtuse triangles
  obtuse_angles = count_obtuse_triangles(cdt);
  std::cout << "\nNumber of obtuse triangles after flips: " << obtuse_angles << std::endl;

  // Draw the triangulation using CGAL's draw function
  CGAL::draw(cdt);

  return 0;
}