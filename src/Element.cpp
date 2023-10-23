#include "Element.hpp"

#include <iostream>

// Public methods

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

std::vector<double> Element::get_x_coordinates() const
{
    std::vector<double> x_coordinates;
    x_coordinates.reserve(n_nodes_);

    for(auto &node : nodes_)
        x_coordinates.push_back(node.get()->get_x());

    return x_coordinates;
}

std::vector<double> Element::get_y_coordinates() const
{
    std::vector<double> y_coordinates;
    y_coordinates.reserve(n_nodes_);

    for(auto &node : nodes_)
        y_coordinates.push_back(node.get()->get_y());

    return y_coordinates;
}

std::vector<double> Element::get_z_coordinates() const
{
    std::vector<double> z_coordinates;
    z_coordinates.reserve(n_nodes_);

    for(auto &node : nodes_)
        z_coordinates.push_back(node.get()->get_z());

    return z_coordinates;
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
