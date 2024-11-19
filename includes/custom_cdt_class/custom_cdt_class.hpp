#ifndef CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H
#define CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Exact_predicates_tag Itag;

#define BOOST_BIND_GLOBAL_PLACEHOLDERS // optional in native ubuntu, removes a warning in wsl

namespace custom_cdt_class {

  template <class Gt, class Tds = CGAL::Default, class Itag = CGAL::Default>

  class Custom_Constrained_Delaunay_triangulation_2

      : public CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag> {

  public:

      using Base = CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>;

      using typename Base::Face_handle;

      using typename Base::Point;

      using typename Base::Vertex_handle;

      using typename Base::Locate_type;

      std::vector<K::FT> steiner_x;
      std::vector<K::FT> steiner_y;



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



      void insert_steiner_x_y(K::FT x, K::FT y) {
        steiner_x.push_back(x);
        steiner_y.push_back(y);
      }

      bool my_is_flippable(const typename Base::Edge& e) {
        Face_handle f1 = e.first; 
        int i = e.second; 
        Face_handle f2 = f1->neighbor(i);
        if (this->is_infinite(f1) || this->is_infinite(f2)) {
            return false;
        }
        else if (this->is_constrained(e)) {
            return false;
        }

        return true;
      }
  };
}


#endif // CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H