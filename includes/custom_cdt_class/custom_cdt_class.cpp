#include "custom_cdt_class.hpp"

// template <class Gt, class Tds, class Itag>
// custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>::Custom_Constrained_Delaunay_triangulation_2(const Gt& gt)

//     : Base(gt) {}

// template <class Gt, class Tds, class Itag>
// custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>::Custom_Constrained_Delaunay_triangulation_2(typename Base::List_constraints& lc, const Gt& gt)

//     : Base(lc, gt) {}

// template <class Gt, class Tds, class Itag>
// template <class InputIterator>
// custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>::Custom_Constrained_Delaunay_triangulation_2(InputIterator it, InputIterator last, const Gt& gt)

//     : Base(it, last, gt) {}

// template <class Gt, class Tds, class Itag>
// typename custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>::Vertex_handle
// custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>::insert_no_flip(const Point& a, Face_handle start) {

//     // Call Ctr::insert without flip_around

//     Vertex_handle va = this->Base::Ctr::insert(a, start); // Directly call Ctr::insert from the base

//     return va;

// }

// template <class Gt, class Tds, class Itag>
// typename custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>::Vertex_handle
// custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>::insert_no_flip(const Point& a, Locate_type lt, Face_handle loc, int li) {

//     Vertex_handle va = this->Base::Ctr::insert(a, lt, loc, li); // Directly call Ctr::insert from the base

//     return va;

// }

// template <class Gt, class Tds, class Itag>
// bool custom_cdt_class::Custom_Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>::my_is_flippable(const typename Base::Edge& e) {
//     Face_handle f1 = e.first; 
//     int i = e.second; 
//     Face_handle f2 = f1->neighbor(i);
//     if (this->is_infinite(f1) || this->is_infinite(f2)) {
//         return false;
//     }
//     else if (this->is_constrained(e)) {
//         return false;
//     }

//     return true;
// }