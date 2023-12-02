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

bool Element::has_thermal_load() const
{
    if (section_.get()->has_temperature())
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

    for (const auto &node : nodes_)
        x_coordinates.push_back(node.get()->get_x());

    return x_coordinates;
}

std::vector<double> Element::get_y_coordinates() const
{
    std::vector<double> y_coordinates;
    y_coordinates.reserve(n_nodes_);

    for (const auto &node : nodes_)
        y_coordinates.push_back(node.get()->get_y());

    return y_coordinates;
}

std::vector<double> Element::get_z_coordinates() const
{
    std::vector<double> z_coordinates;
    z_coordinates.reserve(n_nodes_);

    for (const auto &node : nodes_)
        z_coordinates.push_back(node.get()->get_z());

    return z_coordinates;
}

std::vector<std::size_t> Element::get_nodes() const
{
    std::vector<std::size_t> tags;
    tags.reserve(n_nodes_);

    for (const auto &node : nodes_)
        tags.push_back(node.get()->node_tag());

    return tags;
}

std::vector<std::size_t> Element::local_dofs() const
{
    // Element's nodes' tags
    std::vector<std::size_t> nodes = get_nodes();

    std::vector<std::size_t> dofs;
    dofs.reserve(dof_per_node_ * n_nodes_);

    for (std::size_t i = 0; i < n_nodes_; ++i)
    {
        for (std::size_t j = 0; j < dof_per_node_; ++j)
        {
            std::size_t dof = n_nodes_ * j + i;
            dofs.push_back(dof);
        }
    }

    return dofs;
}

std::vector<std::size_t> Element::local_dofs_by_dof(std::size_t dof) const
{
    std::vector<std::size_t> dofs;
    dofs.reserve(n_nodes_);

    for (std::size_t i = 0; i < n_nodes_; ++i)
        dofs.push_back(i * dof_per_node_ + dof);

    return dofs;
}

std::vector<double>
Element::get_local_displacement_by_dof(const std::vector<double> &q_e,
                                       std::size_t dof) const
{
    std::vector<double> dof_e;
    dof_e.reserve(n_nodes_);

    for (const auto &index : local_dofs_by_dof(dof))
        dof_e.push_back(q_e.at(index));

    return dof_e;
}

std::vector<std::size_t> Element::global_dofs() const
{
    std::vector<std::size_t> nodes = get_nodes();

    std::vector<std::size_t> dofs;
    dofs.reserve(dof_per_node_ * n_nodes_);

    for (std::size_t i = 0; i < n_nodes_; ++i)
    {
        for (std::size_t j = 0; j < dof_per_node_; ++j)
        {
            std::size_t dof = (nodes.at(i) - 1) * dof_per_node_ + j;
            dofs.push_back(dof);
        }
    }

    return dofs;
}

std::vector<std::size_t> Element::global_dofs_by_dof(const std::size_t dof) const
{
    std::vector<std::size_t> global_dofs_vec = global_dofs();

    std::vector<std::size_t> dofs(n_nodes_);

    for (std::size_t i = 0; i < n_nodes_; ++i)
    {
        std::size_t index = dof_per_node_ * i + dof;
        dofs[i] = global_dofs_vec.at(index);
    }

    return dofs;
}

std::vector<double> Element::get_displacements(
    const std::vector<double> &displacements) const
{
    std::vector<double> element_displacement;
    element_displacement.reserve(dof_per_node_ * n_nodes_);

    for (auto const &index : global_dofs())
        element_displacement.push_back(displacements.at(index));

    return element_displacement;
}

std::vector<double> Element::get_w_displacements(
    const std::vector<double> &displacements) const
{
    std::vector<double> w_displacement;
    w_displacement.reserve(n_nodes_);

    for (auto const &index : global_dofs_by_dof(2))
        w_displacement.push_back(displacements.at(index));

    return w_displacement;
}

// Check if element has any node with boundary condition
bool Element::has_boundary_node() const
{
    bool boundary = false;

    for (auto const &node : nodes_)
    {
        if (node.get()->has_boundary())
        {
            boundary = true;
            return boundary;
        }
    }

    return boundary;
}

std::vector<std::size_t> Element::boundary_dofs() const
{
    std::vector<std::size_t> boundary_dofs{};

    if (has_boundary_node())
    {
        for (auto const &node : nodes_)
        {
            if (node.get()->has_boundary())
            {
                std::size_t node_tag = node.get()->node_tag();

                BoundaryPtr boundary = node.get()->get_boundary();
                for (std::size_t i = 0; i < dof_per_node_; ++i)
                {
                    bool has_constraint = boundary.get()->imposed_dofs().at(i);

                    if (has_constraint)
                    {
                        std::size_t dof = (node_tag - 1) * dof_per_node_ + i;
                        boundary_dofs.push_back(dof);
                    }
                }
            }
        }
    }

    return boundary_dofs;
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