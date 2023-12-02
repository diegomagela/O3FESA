#include "Triplet.hpp"

std::ostream &operator<<(std::ostream &os, const Triplet &triplet)
{
    os << triplet.row_ << ','
       << triplet.col_ << ','
       << triplet.value_;

    return os;
}