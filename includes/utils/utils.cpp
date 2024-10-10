#include "utils.hpp"

string utils::simple_or_not(const Polygon_2 & p)
{
    if (p.is_simple())
        return "Simple Polygon";
    else 
        return "Not Simple Polygon";
}


