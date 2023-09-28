#ifndef S9_H
#define S9_H

#include <iostream>

#include "Element.hpp"
#include "Node.hpp"

// First shear deformation theory second-order quadrangle lagrange element
class S9 : public Element
{
public:
    S9(const std::size_t tag,
       const std::vector<NodePtr> nodes) : Element(type_,
                                                   tag,
                                                   n_nodes_,
                                                   dof_per_node_,
                                                   nodes){};
    // stiffness_matrix () override
    // mass_matrix () override

    // The rule of five
    S9() = default;
    S9(S9 const &) = default;
    S9 &operator=(S9 const &) = default;
    S9(S9 &&) = default;
    S9 &operator=(S9 &&) = default;

private:
    inline static std::string type_{"S9"};
    inline static constexpr std::size_t n_nodes_ = 9;
    inline static constexpr std::size_t dof_per_node_ = 5;
};

#endif // S9_H