#ifndef ELEMENT_H
#define ELEMENT_H

#include <cstddef>
#include <memory>
#include <vector>

#include "DLoad.hpp"
#include "Node.hpp"
#include "Section.hpp"

typedef std::shared_ptr<Node> NodePtr;
typedef std::shared_ptr<Section> SectionPtr;

/*
    TODO

    For better performance (?), consider storing only a vector of nodes' tags,
    instead of a vector of Nodes.

*/

// Base element class for all shell other elements
class Element
{
public:
    explicit Element(const std::string type,
                     const std::size_t tag,
                     const std::size_t n_nodes,
                     const std::size_t dof_per_node,
                     const std::vector<NodePtr> nodes) : type_(type),
                                                         tag_(tag),
                                                         n_nodes_(n_nodes),
                                                         dof_per_node_(dof_per_node),
                                                         nodes_(nodes){};

    // The rule of five
    Element() = default;
    Element(Element const &) = default;
    Element &operator=(Element const &) = default;
    Element(Element &&) = default;
    Element &operator=(Element &&) = default;

    void print();

    // Setter
    inline void set_dload(DLoad dload) { dload_ = dload; };
    inline void set_section(SectionPtr section) { section_ = section; };

    bool has_load() { return dload_.has_load(); };

    // Pure virtual functions to be implemented

    // stiffness_matrix() = 0;
    // mass_matrix() = 0;

    friend std::ostream &operator<<(std::ostream &os, const Element &element);

private:
    std::string type_{};           // Element type
    std::size_t tag_{};            // Element's tag
    std::size_t n_nodes_{};        // Number of nodes
    std::size_t dof_per_node_{};   // Number of degree of freedom per node
    std::vector<NodePtr> nodes_{}; // List of element's nodes
    DLoad dload_{};                // Area distributed loading
    SectionPtr section_{};
};

#endif // ELEMENT_H
