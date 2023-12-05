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

std::vector<std::size_t> Element::get_nodes_tags() const
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
    std::vector<std::size_t> nodes = get_nodes_tags();

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
    std::vector<std::size_t> nodes = get_nodes_tags();

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

std::vector<double> Element::get_displacement() const
{
    std::vector<double> displacements;
    displacements.reserve(total_dof());

    for (auto const &node : get_nodes())
    {
        std::vector<double> node_displacement = node.get()->get_displacements();
        vct::append_vector(displacements, node_displacement);
    }

    return displacements;
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

std::vector<Triplet>
Element::matrix_to_triplet(const Eigen::MatrixXd &matrix) const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> matrix_triplet;
    matrix_triplet.reserve(estimation_nnz);

    std::vector<std::size_t> element_global_dofs = global_dofs();
    std::vector<std::size_t> element_local_dofs = local_dofs();
    std::vector<std::size_t> bc_dofs = boundary_dofs();

    for (std::size_t i = 0; i < total_dof(); ++i)
    {
        for (std::size_t j = 0; j < total_dof(); ++j)
        {
            // Local indexes
            std::size_t row_local = element_local_dofs.at(i);
            std::size_t col_local = element_local_dofs.at(j);

            // Global indexes
            std::size_t row_global = element_global_dofs.at(i);
            std::size_t col_global = element_global_dofs.at(j);

            double value = matrix(row_local, col_local);

            if (check_boundary_dof_matrix(bc_dofs, row_global, col_global))
            {
                if (row_global == col_global)
                    value = matrix.mean();

                else
                    value = 0;
            }

            if (std::abs(value) > 0)
            {
                Triplet t(row_global, col_global, value);
                matrix_triplet.push_back(t);
            }
        }
    }

    return matrix_triplet;
}

std::vector<Triplet>
Element::vector_to_triplet(const Eigen::VectorXd &vector) const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> load;
    load.reserve(estimation_nnz);

    std::vector<std::size_t> element_global_dofs = global_dofs();
    std::vector<std::size_t> element_local_dofs = local_dofs();
    std::vector<std::size_t> bc_dofs = boundary_dofs();

    for (std::size_t i = 0; i < total_dof(); ++i)
    {
        // Local index
        std::size_t row_local = element_local_dofs.at(i);

        // Global index
        std::size_t row_global = element_global_dofs.at(i);

        double value = vector(row_local);

        if (check_boundary_dof_vector(bc_dofs, row_global))
            value = 0;

        if (std::abs(value) > 0)
        {
            Triplet t(row_global, 0, value);
            load.push_back(t);
        }
    }

    return load;
}

bool Element::check_boundary_dof_vector(const std::vector<std::size_t> &dofs,
                                        const std::size_t row) const
{
    bool check = false;

    if (vct::vector_has_element(dofs, row))
        check = true;

    return check;
}

bool Element::check_boundary_dof_matrix(const std::vector<std::size_t> &dofs,
                                        const std::size_t row,
                                        const std::size_t col) const
{
    bool check = false;

    if (vct::vector_has_element(dofs, row) or vct::vector_has_element(dofs, col))
        check = true;

    return check;
}

std::vector<Triplet> Element::element_load_triplet() const
{
    const std::size_t size = pressure_load_triplet().size() +
                             thermal_load_vector().size();

    std::vector<Triplet> element_load;
    element_load.reserve(size);

    vct::append_vector(element_load, pressure_load_triplet());
    vct::append_vector(element_load, thermal_load_triplet());

    return element_load;
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