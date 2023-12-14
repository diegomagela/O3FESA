#include "Element.hpp"

#include <iostream>

// Public methods

bool Element::has_load() const
{
    bool load{false};

    if (dload_)
        load = true;

    return load;
}

bool Element::has_thermal_load() const
{
    bool thermal_load{false};

    if (section_.get()->has_temperature())
        thermal_load = true;

    return thermal_load;
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

std::vector<double> Element::get_u_displacement() const
{
    std::vector<double> displacements;
    displacements.reserve(dof_per_node_);

    for (auto const &node : get_nodes())
    {
        double node_u_displacement = node.get()->get_u_displacements();
        displacements.push_back(node_u_displacement);
    }

    return displacements;
}

std::vector<double> Element::get_v_displacement() const
{
    std::vector<double> displacements;
    displacements.reserve(dof_per_node_);

    for (auto const &node : get_nodes())
    {
        double node_v_displacement = node.get()->get_v_displacements();
        displacements.push_back(node_v_displacement);
    }

    return displacements;
}

std::vector<double> Element::get_w_displacement() const
{
    std::vector<double> displacements;
    displacements.reserve(dof_per_node_);

    for (auto const &node : get_nodes())
    {
        double node_w_displacement = node.get()->get_w_displacements();
        displacements.push_back(node_w_displacement);
    }

    return displacements;
}

std::vector<double> Element::get_phi_x_displacement() const
{
    std::vector<double> displacements;
    displacements.reserve(dof_per_node_);

    for (auto const &node : get_nodes())
    {
        double node_phi_x_displacement = node.get()->get_phi_x_displacements();
        displacements.push_back(node_phi_x_displacement);
    }

    return displacements;
}

std::vector<double> Element::get_phi_y_displacement() const
{
    std::vector<double> displacements;
    displacements.reserve(dof_per_node_);

    for (auto const &node : get_nodes())
    {
        double node_phi_y_displacement = node.get()->get_phi_y_displacements();
        displacements.push_back(node_phi_y_displacement);
    }

    return displacements;
}

std::vector<double> Element::get_displacements() const
{
    std::vector<double> displacements{};
    displacements.reserve(total_dof());

    std::vector<double> u_e = get_u_displacement();
    std::vector<double> v_e = get_v_displacement();
    std::vector<double> w_e = get_w_displacement();
    std::vector<double> phi_x_e = get_phi_x_displacement();
    std::vector<double> phi_y_e = get_phi_y_displacement();

    // >>>
    displacements.insert(displacements.end(),
                         u_e.begin(),
                         u_e.end());

    displacements.insert(displacements.end(),
                         v_e.begin(),
                         v_e.end());

    displacements.insert(displacements.end(),
                         w_e.begin(),
                         w_e.end());

    displacements.insert(displacements.end(),
                         phi_x_e.begin(),
                         phi_x_e.end());

    displacements.insert(displacements.end(),
                         phi_y_e.begin(),
                         phi_y_e.end());

    return displacements;
}

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

std::vector<std::size_t> Element::local_boundary_dofs() const
{
    std::vector<std::size_t> boundary_dofs{};

    if (has_boundary_node())
    {
        for (std::size_t i = 0; i < n_nodes_; ++i)
        {
            NodePtr node = nodes_.at(i);

            if (node.get()->has_boundary())
            {
                BoundaryPtr boundary = node.get()->get_boundary();

                for (std::size_t j = 0; j < dof_per_node_; ++j)
                {
                    bool has_constraint = boundary.get()->imposed_dofs().at(j);

                    if (has_constraint)
                    {
                        std::size_t dof = i * dof_per_node_ + j;
                        boundary_dofs.push_back(dof);
                    }
                }
            }
        }
    }

    return boundary_dofs;
}

std::vector<std::size_t> Element::global_boundary_dofs() const
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
    std::vector<std::size_t> bc_dofs = global_boundary_dofs();

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

            if (check_boundary_dof_matrix(row_global, col_global))
            {
                if (row_global == col_global)
                    value = 1;
                // value = matrix.mean();

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
    std::vector<std::size_t> bc_dofs = global_boundary_dofs();

    for (std::size_t i = 0; i < total_dof(); ++i)
    {
        // Local index
        std::size_t row_local = element_local_dofs.at(i);

        // Global index
        std::size_t row_global = element_global_dofs.at(i);

        double value = vector(row_local);

        if (check_boundary_dof_vector(row_global))
            value = 0;

        if (std::abs(value) > 0)
        {
            Triplet t(row_global, 0, value);
            load.push_back(t);
        }
    }

    return load;
}

bool Element::check_boundary_dof_vector(const std::size_t row) const
{
    bool check = false;

    if (vct::vector_has_element(global_boundary_dofs(), row))
        check = true;

    return check;
}

bool Element::check_boundary_dof_matrix(const std::size_t row,
                                        const std::size_t col) const
{
    bool check = false;

    if (vct::vector_has_element(global_boundary_dofs(), row) or
        vct::vector_has_element(global_boundary_dofs(), col))
        check = true;

    return check;
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