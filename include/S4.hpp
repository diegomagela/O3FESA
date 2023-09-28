#ifndef S4_H
#define S4_H

#include <iostream>

#include "Element.hpp"
#include "Node.hpp"

// First-order quadrangle lagrange element
class S4 : public Element
{
public:
    S4(const std::size_t tag,
       const std::vector<NodePtr> nodes) : Element(type_,
                                                   tag,
                                                   n_nodes_,
                                                   dof_per_node_,
                                                   nodes){};

private:
    inline static std::string type_{"S4"};
    inline static constexpr std::size_t n_nodes_ = 4;
    inline static constexpr std::size_t dof_per_node_ = 5;
};

#endif // S4_H