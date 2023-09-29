#include "Element.hpp"

#include <iostream>

void Element::print() const
{
    std::cout << *this << std::endl;
}

std::ostream &operator<<(std::ostream &os, const Element &element)
{
    os << element.type_ << '\t'
       << element.tag_ << '\t';

    for (auto const &node : element.nodes_)
        os << node.get()->node_tag() << '\t';

    return os;
}
