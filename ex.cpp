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

//
class obt_point {
  public:
    int obt_count;
    Point insrt_pt;

    obt_point(int count, Point pt) {
      obt_count = count;
      insrt_pt = pt;
    }
};
//

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

    if (!(v1.x() == p1_point.x() && v1.y() == p1_point.y() && v2.x() == p2_point.x() && v2.y() == p2_point.y())) {
      continue;
    }
    else 
      std::cout << "Found the edge: (" << p1->point() << ") - (" << p2->point() << ")" << std::endl;

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

obt_point insert_circumcenter(CDT& cdt, CDT::Face_handle f1) {

  // Calculate the circumcenter of the triangle
  Point a = f1->vertex(0)->point();
  Point b = f1->vertex(1)->point();
  Point c = f1->vertex(2)->point();
  Point pericenter = CGAL::circumcenter(a, b, c);

  // Check if the inserted vertex is inside the convex hull
  CDT::Face_handle located_face = cdt.locate(pericenter);
  if (!cdt.is_infinite(located_face)) {
    cdt.insert_no_flip(pericenter);
  }

  // return count_obtuse_triangles(copy);
  obt_point ret(count_obtuse_triangles(cdt), pericenter);
  return ret;
}

obt_point insert_centroid(CDT& cdt, CDT::Face_handle f1) {

  // Calculate the centroid of the triangle
  Point a = f1->vertex(0)->point();
  Point b = f1->vertex(1)->point();
  Point c = f1->vertex(2)->point();
  Point centroid = CGAL::centroid(a, b, c);

  // Insert the centroid
  cdt.insert_no_flip(centroid);


  // count_obtuse_triangles(copy);
  obt_point ret(count_obtuse_triangles(cdt), centroid);
  return ret;
}

int merge_obtuse(CDT& cdt, CDT::Face_handle f1) {
  CDT copy(cdt);

  // Get the vertices of the triangle
  Point p1 = f1->vertex(0)->point();
  Point p2 = f1->vertex(1)->point();
  Point p3 = f1->vertex(2)->point();

  // Check if the neighbors of the triangle are obtuse and merge the obtuse triangles
  CDT::Face_handle neigh1 = f1->neighbor(0);
  CDT::Face_handle neigh2 = f1->neighbor(1);
  CDT::Face_handle neigh3 = f1->neighbor(2);

  // Check if the neighbors are obtuse
  bool obt1 = has_obtuse_angle(neigh1);

  if (has_obtuse_angle(neigh1)) {
    CGAL::draw(copy);
    int i = f1->index(neigh1);
    

    // DEN EXO IDEA POS NA KANO MERGE TA TRIGONA lol


    // Get the point that is exactly in the center of the shared edge
    CDT::Vertex_handle v1 = f1->vertex((i + 1) % 3);
    CDT::Vertex_handle v2 = f1->vertex((i + 2) % 3);

    // Calculate the midpoint of the shared edge
    Point midpoint = CGAL::midpoint(v1->point(), v2->point());
    std::cout << "Midpoint of the shared edge: " << midpoint << std::endl;
    copy.insert_no_flip(midpoint);
    CGAL::draw(copy);
  }



  return count_obtuse_triangles(copy);
}


obt_point insert_mid(CDT& cdt, CDT::Face_handle f1) {
  CDT copy(cdt);

  // Calculate the circumcenter of the triangle
  Point a = f1->vertex(0)->point();
  Point b = f1->vertex(1)->point();
  Point c = f1->vertex(2)->point();
  Point mid;

  K::FT l0 = CGAL::squared_distance(a, b);
  K::FT l1 = CGAL::squared_distance(a, c);
  K::FT l2 = CGAL::squared_distance(b, c);

  // Calculate the midpoint of edge with the longest length
  Point midpoint;
  if (l0 >= l1 && l0 >= l2) {
    midpoint = CGAL::midpoint(a, b);
  }
  else if (l1 >= l2) {
    midpoint = CGAL::midpoint(a, c);
  }
  else {
    midpoint = CGAL::midpoint(b, c);
  }

  // std::cout << "midpoint: " << midpoint << " | points: a:" << a << " b: " << b << std::endl;
  cdt.insert_no_flip(midpoint);
  int obt = count_obtuse_triangles(cdt);
  obt_point ret = obt_point(obt, midpoint);
  return ret;
}

void steiner_insertion(CDT& cdt) {
  int init_obtuse_count = count_obtuse_triangles(cdt);
  std::cout << "Initial obtuse count: " << init_obtuse_count << std::endl;
  Point a;
  obt_point best_steiner(-1, a);

  // Iterate the faces of the cdt
  for (CDT::Finite_faces_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++) {

    if (has_obtuse_angle(f)) {
      
      CDT copy(cdt);
      // Insert the circumcenter if possible
      obt_point calc_insert_mid = insert_mid(copy, f);
      if (best_steiner.obt_count == -1) {
        best_steiner = calc_insert_mid;
      }
      else if (best_steiner.obt_count > calc_insert_mid.obt_count) {
        best_steiner = calc_insert_mid;
      }
      // std::cout << "Obtuse triangles after inserting the edge midpoint: " << count_obtuse_triangles(copy) << std::endl;

      CDT copy1(cdt);
      // Insert the circumcenter if possible
      obt_point calc_insert_centr = insert_centroid(copy1, f);
      if (best_steiner.obt_count > calc_insert_centr.obt_count) {
        best_steiner = calc_insert_centr;
      }
      // std::cout << "Obtuse triangles after inserting the centroid: " << count_obtuse_triangles(copy1) << std::endl;

      CDT copy2(cdt);
      // Insert the circumcenter if possible
      obt_point calc_insert_circ = insert_circumcenter(copy2, f);
      if (best_steiner.obt_count > calc_insert_circ.obt_count) {
        best_steiner = calc_insert_circ;
      }
      // std::cout << "Obtuse triangles after inserting the circumcenter: " << count_obtuse_triangles(copy2) << std::endl;
    
    }
  }
  // if (best_steiner.obt_count <= count_obtuse_triangles(cdt)) {
  // CGAL::draw(cdt);
  std::cout << "largest edge midpoint: " << best_steiner.insrt_pt << std::endl;
  cdt.insert_no_flip(best_steiner.insrt_pt);
  // CGAL::draw(cdt);
  // }

  std::cout << "Final obtuse count: " << count_obtuse_triangles(cdt) << std::endl;
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
  std::cout << "Number of obtuse triangles before the flips: " << count_obtuse_triangles(cdt) << std::endl;
  
  // Make flips
  make_flips(cdt);

  // Insert Steiner points
  std::cout << "\n\n\nSteiner points insertion:\n";
  for (int i = 0 ; i < 200 ; i++) {
    steiner_insertion(cdt);
  }

  // Count the obtuse triangles
  std::cout << "Number of obtuse triangles after the flips: " << count_obtuse_triangles(cdt) << std::endl;

  // Draw the triangulation using CGAL's draw function
  CGAL::draw(cdt);

  return 0;
}