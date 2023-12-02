#include "CLoad.hpp"

std::vector<std::size_t> CLoad::local_dofs() const
{
    std::vector<std::size_t> dofs;
    dofs.reserve(6); // Maximum number of dofs possible in each node

    for (std::size_t i = 0; i < loading_dofs_.size(); ++i)
    {
        if (loading_dofs_[i])
            dofs.push_back(i);
    }

    return dofs;
}

std::vector<std::size_t> CLoad::global_dofs() const
{
    // If 6 degrees of freedom per node are considered?
    constexpr std::size_t n_dof = 5;

    std::vector<std::size_t> dofs;
    dofs.reserve(6); // Maximum number of dofs possible in each node

    for (std::size_t i = 0; i < loading_dofs_.size(); ++i)
    {
        if (loading_dofs_[i])
        {
            std::size_t dof = (node_tag() - 1) * n_dof + i;
            dofs.push_back(dof);
        }
    }

    return dofs;
}

std::vector<double> CLoad::load_values() const
{
    std::vector<double> loads;
    loads.reserve(6); // Maximum number of dofs possible in each node

    for (std::size_t i = 0; i < loading_dofs_.size(); ++i)
    {
        if (loading_dofs_[i])
            loads.push_back(loading_values_.at(i));
    }

    return loads;
}