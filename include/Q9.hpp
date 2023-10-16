#ifndef Q9_H
#define Q9_H

#include <iostream>

#include "Element.hpp"
#include "Node.hpp"

// First shear deformation theory second-order quadrangle lagrange element
class Q9 : public Element
{
public:
    Q9(const std::size_t tag) : Element(type_,
                                        tag,
                                        n_nodes_,
                                        dof_per_node_){};

    // The rule of five
    Q9() = default;
    Q9(Q9 const &) = default;
    Q9 &operator=(Q9 const &) = default;
    Q9(Q9 &&) = default;
    Q9 &operator=(Q9 &&) = default;

    // stiffness_matrix () override
    // mass_matrix () override

private:
    std::vector<double> shape_functions(const double xi, const double eta) const;
    std::vector<double> shape_functions_d_xi(const double xi,
                                             const double eta) const;
    std::vector<double> shape_functions_d_eta(const double xi,
                                              const double eta) const;

private:
    inline static std::string type_{"Q9"};
    inline static constexpr std::size_t n_nodes_ = 9;
    inline static constexpr std::size_t dof_per_node_ = 5;
};

#endif // Q9_H