#ifndef UTILS_HPP
#define UTILS_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon_2;

using std::string;

namespace utils {
    string simple_or_not(const Polygon_2 & p);
}

#endif
