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

    // Selectors

    // Returns if has any boundary dof
    inline bool has_boundary() const
    {
        bool bc = true;

        if (imposed_dofs_.empty())
            bc = false;

        return bc;
    }

    // Returns imposed dofs vector
    inline std::vector<bool> imposed_dofs() const { return imposed_dofs_; }

    // Modifiers

    // Set a new imposed dof 
    inline void set_imposed_dof(const std::size_t dof)
    {
        imposed_dofs_.at(dof) = true;
    }

private:
    std::size_t node_tag_{};
    std::vector<bool> imposed_dofs_{};
    std::vector<double> imposed_values_{};
};

#endif // BOUNDARY_HPP
