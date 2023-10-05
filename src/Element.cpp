#include "Element.hpp"

#include <iostream>

bool Element::has_load() const
{
    if (dload_)
        return true;

    else
        return false;
}

void Element::print() const
{
    std::cout << *this << std::endl;
}

std::ostream &operator<<(std::ostream &os, const Element &element)
{
    os << element.tag_ << '\t'
       << element.type_ << '\t'
       << element.has_load() << '\t';


    for (auto const &node : element.nodes_)
        os << node.get()->node_tag() << '\t';

    return os;
}
