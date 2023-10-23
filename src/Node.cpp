#include "Node.hpp"

#include <iostream>

bool Node::has_boundary() const
{
    if (boundary_)
        return true;

    else
        return false;
}

bool Node::has_cload() const
{
    if (cload_)
        return true;

    else
        return false;
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