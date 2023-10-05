#ifndef ELEMENT_H
#define ELEMENT_H

#include <cstddef>
#include <memory>
#include <vector>

#include "DLoad.hpp"
#include "Node.hpp"
#include "Section.hpp"

/*
    TODO

    For better performance (?), consider storing only a vector of nodes' tags,
    instead of a vector of Nodes.

*/

// Base element class for all shell other elements
class Element
{
private:
    typedef std::shared_ptr<Node> NodePtr;
    typedef std::shared_ptr<Section> SectionPtr;
    typedef std::shared_ptr<DLoad> DLoadPtr;

public:
    explicit Element(const std::string type,
                     const std::size_t tag,
                     const std::size_t n_nodes,
                     const std::size_t dof_per_node) : type_(type),
                                                       tag_(tag),
                                                       n_nodes_(n_nodes),
                                                       dof_per_node_(dof_per_node){};

    // The rule of five
    Element() = default;
    Element(Element const &) = default;
    Element &operator=(Element const &) = default;
    Element(Element &&) = default;
    Element &operator=(Element &&) = default;

    // Selectors
    bool has_load() const;
    void print() const;

    // Modifiers
    inline void set_nodes(std::vector<NodePtr> nodes) { nodes_ = nodes; }
    inline void set_dload(DLoadPtr dload) { dload_ = dload; }
    inline void set_section(SectionPtr section) { section_ = section; }

    // Pure virtual functions to be implemented

    // stiffness_matrix() = 0;
    // mass_matrix() = 0;

    // Friend
    friend std::ostream &operator<<(std::ostream &os, const Element &element);

private:
    std::string type_{};           // Element type
    std::size_t tag_{};            // Element's tag
    std::size_t n_nodes_{};        // Number of nodes
    std::size_t dof_per_node_{};   // Number of degree of freedom per node
    std::vector<NodePtr> nodes_{}; // List of element's nodes
    DLoadPtr dload_;               // Area distributed loading
    SectionPtr section_;
};

#endif // ELEMENT_H
