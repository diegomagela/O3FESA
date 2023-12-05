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
    inline bool has_boundary() const
    {
        bool bc = true;

        if (imposed_dofs_.empty())
            bc = false;

        return bc;
    }

    inline std::vector<bool> imposed_dofs() const { return imposed_dofs_; }

private:
    std::size_t node_tag_{};
    std::vector<bool> imposed_dofs_{};
    std::vector<double> imposed_values_{};
};

#endif // BOUNDARY_HPP
