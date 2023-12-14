#ifndef ELEMENT_H
#define ELEMENT_H

#include <cstddef>
#include <memory>
#include <vector>

#include <Eigen/Dense>

#include "DLoad.hpp"
#include "Node.hpp"
#include "Section.hpp"
#include "VectorUtil.hpp"
#include "Triplet.hpp"

/*
    TODO

    Is it better consider storing only a vector of nodes' tags, instead of a
    vector of Nodes? I do not think so, but I must benchmark this.

*/

// Base element class for all other plate/shell elements
class Element
{
private:
    typedef std::shared_ptr<Node> NodePtr;
    typedef std::shared_ptr<Section> SectionPtr;
    typedef std::shared_ptr<DLoad> DLoadPtr;
    typedef std::shared_ptr<Boundary> BoundaryPtr;

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

    //// Modifiers

    inline void set_nodes(std::vector<NodePtr> &nodes) { nodes_ = nodes; }
    inline void set_dload(DLoadPtr dload) { dload_ = dload; }
    inline void set_section(SectionPtr section) { section_ = section; }

    //// Selectors

    // Inline functions

    // Total number of dofs of the element
    inline std::size_t total_dof() const { return dof_per_node_ * n_nodes_; }
    // Returns the pressure load value
    inline double load_value() const { return dload_.get()->load_value(); }
    // Return elements' nodes objects
    inline std::vector<NodePtr> get_nodes() const { return nodes_; }
    // Returns the element's section
    inline SectionPtr get_section() const { return section_; }

    // Non-inline functions

    // Returns if the element has a pressure load
    bool has_load() const;
    // Returns if the element has a temperature gradient
    bool has_thermal_load() const;
    // Return if any node of the element has a boundary condition
    bool has_boundary_node() const;
    // Returns the x coordinates of all nodes in the element
    std::vector<double> get_x_coordinates() const;
    // Returns the y coordinates of all nodes in the element
    std::vector<double> get_y_coordinates() const;
    // Returns the z coordinates of all nodes in the element
    std::vector<double> get_z_coordinates() const;
    // Return element's w displacement vector for all nodes
    std::vector<double> get_u_displacement() const;
    // Return element's w displacement vector for all nodes
    std::vector<double> get_v_displacement() const;
    // Return element's w displacement vector for all nodes
    std::vector<double> get_w_displacement() const;
    // Return element's w displacement vector for all nodes
    std::vector<double> get_phi_x_displacement() const;
    // Return element's w displacement vector for all nodes
    std::vector<double> get_phi_y_displacement() const;
    // Return element's displacement vector
    std::vector<double> get_displacements() const;
    // Return the global indexes of the boundary conditions in the element
    std::vector<std::size_t> global_boundary_dofs() const;
    // Return the local indexes of the boundary conditions in the element
    std::vector<std::size_t> local_boundary_dofs() const;
    // Returns the tags of all nodes
    std::vector<std::size_t> get_nodes_tags() const;
    // Returns the element local indexes
    std::vector<std::size_t> local_dofs() const;
    // Returns the element global indexes
    std::vector<std::size_t> global_dofs() const;
    // Convert dense matrix to triplet format applying boundary conditions
    std::vector<Triplet> matrix_to_triplet(const Eigen::MatrixXd &matrix) const;
    // Convert dense vector to triplet format applying boundary conditions
    std::vector<Triplet> vector_to_triplet(const Eigen::VectorXd &vector) const;
    // Check if dof in a vector has boundary condition
    bool check_boundary_dof_vector(const std::size_t row) const;
    // Check if dof in a matrix has boundary condition
    bool check_boundary_dof_matrix(const std::size_t row,
                                   const std::size_t col) const;

    // Virtual functions

    // Triplet linear stiffness matrix
    virtual std::vector<Triplet> linear_stiffness_triplet() const = 0;
    // Triplet linear stiffness matrix
    virtual std::vector<Triplet> nonlinear_stiffness_triplet() const = 0;
    // Triplet external load vector
    virtual std::vector<Triplet> external_load_triplet() const = 0;
    // Triplet internal load vector: K(q)q
    virtual std::vector<Triplet> internal_load_triplet() const = 0;
    // Triplet tangent stiffness matrix
    virtual std::vector<Triplet> tangent_stiffness_triplet() const = 0;

    // Friend functions

    // Friend operator <<
    friend std::ostream &operator<<(std::ostream &os, const Element &element);

private:
    std::string type_{};           // Element type
    std::size_t tag_{};            // Element tag
    std::size_t n_nodes_{};        // Number of nodes
    std::size_t dof_per_node_{};   // Number of degree of freedom per node
    std::vector<NodePtr> nodes_{}; // List of element nodes
    DLoadPtr dload_;               // Area distributed loading
    SectionPtr section_;           // Element section
};

#endif // ELEMENT_H
