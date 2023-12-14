#include "Node.hpp"

#include <iostream>

bool Node::has_boundary() const
{
    bool bc{false};

    if (boundary_)
        bc = true;

    return bc;
}

bool Node::has_cload() const
{
    bool load{false};

    if (cload_)
        load = true;

    return load;
}

std::ostream &operator<<(std::ostream &os, const Node &node)
{
    os << node.tag_ << '\t'
       << node.x_ << '\t'
       << node.y_ << '\t'
       << node.y_ << '\t'
       << node.has_boundary() << '\t'
       << node.has_cload() << '\t';

    return os;
}

void Node::print() const
{
    std::cout << *this << std::endl;
}