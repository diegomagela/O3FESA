#ifndef Q4_H
#define Q4_H

#include <iostream>

#include "Element.hpp"
#include "Node.hpp"

// First-order quadrangle lagrange element
class Q4 : public Element
{
public:
    Q4(const std::size_t tag) : Element(type_,
                                        tag,
                                        n_nodes_,
                                        dof_per_node_){};

private:
    inline static std::string type_{"Q4"};
    inline static constexpr std::size_t n_nodes_ = 4;
    inline static constexpr std::size_t dof_per_node_ = 5;
};

#endif // Q4_H