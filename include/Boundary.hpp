#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <vector>

class Boundary
{
public:
    Boundary(std::size_t node_tag,
             std::vector<bool> imposed_dofs,
             std::vector<double> imposed_values) : node_tag_(node_tag),
                                                   imposed_dofs_(imposed_dofs),
                                                   imposed_values_(imposed_values){};

    // The rule of five
    Boundary() = default;
    Boundary(Boundary const &) = default;
    Boundary &operator=(Boundary const &) = default;
    Boundary(Boundary &&) = default;
    Boundary &operator=(Boundary &&) = default;

    bool has_boundary()
    {
        if(imposed_dofs_.empty())
            return false;

        else
            return true;
    }

private:
    std::size_t node_tag_{};
    std::vector<bool> imposed_dofs_{};
    std::vector<double> imposed_values_{};
};

#endif // BOUNDARY_HPP
