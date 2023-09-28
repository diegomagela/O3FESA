#include "Node.hpp"

#include <iostream>

std::ostream &operator<<(std::ostream &os, const Node &node)
{
    os << node.tag_ << '\t'
       << node.x_ << '\t'
       << node.y_ << '\t'
       << node.y_ << '\t';

    return os;
}

void Node::print() const
{
    std::cout << *this << std::endl;
}