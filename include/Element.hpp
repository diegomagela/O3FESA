#ifndef ELEMENT_H
#define ELEMENT_H

#include <cstddef>
#include <memory>
#include <vector>

#include <Eigen/Dense>

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

    //// Selectors

    // Total number of dofs of the element
    inline std::size_t total_dof() const { return dof_per_node_ * n_nodes_; }

    // Returns if the element has a pressure load
    bool has_load() const;

    bool has_thermal_load() const;

    // Returns the pressure load value
    inline double load_value() const { return dload_.get()->load_value(); }

    // Returns the x coordinates of all nodes in the element
    std::vector<double> get_x_coordinates() const;

    // Returns the y coordinates of all nodes in the element
    std::vector<double> get_y_coordinates() const;

    // Returns the z coordinates of all nodes in the element
    std::vector<double> get_z_coordinates() const;

    // Returns the tags of all element' nodes
    std::vector<std::size_t> get_nodes() const;

    // Returns the element local indexes
    std::vector<std::size_t> local_dofs() const;

    // Returns the element local indexes by dof
    // u:       1
    // v:       2
    // w:       3
    // phi_x:   4
    // phi_y:   5
    std::vector<std::size_t> local_dofs_by_dof(const std::size_t dof) const;

    inline std::vector<double>
    get_local_w_displacement(const std::vector<double> &q_e) const
    {
        return get_local_displacement_by_dof(q_e, 2);
    }

    // Returns the element global indexes
    std::vector<std::size_t> global_dofs() const;

    // Return element displacement vector for all dofs
    // displacements: total displacements
    std::vector<double>
    get_displacements(const std::vector<double> &displacements) const;

    // Return the element's w displacement vector
    // displacements: total displacements
    std::vector<double>
    get_w_displacements(const std::vector<double> &displacements) const;

    // Return if the any node of the element has a boundary condition
    bool has_boundary_node() const;

    // Return the global indexes of the boundary conditions in the element
    std::vector<std::size_t> boundary_dofs() const;

    // Returns the element section
    inline SectionPtr get_section() const { return section_; }

    // Prints to the console the element information
    void print() const;

    // Modifiers

    inline void set_nodes(std::vector<NodePtr> &nodes) { nodes_ = nodes; }
    inline void set_dload(DLoadPtr dload) { dload_ = dload; }
    inline void set_section(SectionPtr section) { section_ = section; }

    //// Pure virtual functions

    /// Jacobian transformation methods

    virtual std::vector<double> jacobian_matrix(const double xi,
                                                const double eta) const = 0;

    virtual double jacobian(const double xi, const double eta) const = 0;

    virtual std::vector<double>
    jacobian_inverse_matrix(const double xi,
                            const double eta) const = 0;

    // Matrices //

    virtual Eigen::MatrixXd shape_functions_matrix(const double xi,
                                                   const double eta) const = 0;

    // Dense stiffness matrix
    virtual Eigen::MatrixXd linear_stiffness_matrix() const = 0;
    virtual Eigen::MatrixXd
    nonlinear_stiffness_matrix(const std::vector<double> &w_e) const = 0;

    //// TESTING
    virtual Eigen::MatrixXd
    tangent_stiffness_matrix(const std::vector<double> &q_e) const = 0;

    // f_int = K^e(d^e)d^e
    // d_e: total displacement within the element
    virtual Eigen::VectorXd
    internal_force(const std::vector<double> &d_e) const = 0;

    // Dense mass matrix
    // virtual Eigen::MatrixXd mass_matrix() const = 0;

    // Vectors //

    // Dense pressure load vector
    virtual Eigen::VectorXd pressure_load_vector() const = 0;

    // Dense thermal load vector
    virtual Eigen::VectorXd thermal_load_vector() const = 0;

    // Friend
    friend std::ostream &operator<<(std::ostream &os, const Element &element);

private:
    // Returns the element global indexes according to dof number
    // u =      0
    // v =      1
    // w =      2
    // phi_x =  3
    // phi_y =  4
    std::vector<std::size_t> global_dofs_by_dof(const std::size_t dof) const;

    // Returns the dof element displacement vector
    // u =      0
    // v =      1
    // w =      2
    // phi_x =  3
    // phi_y =  4
    std::vector<double>
    get_local_displacement_by_dof(const std::vector<double> &q_e,
                                  const std::size_t dof) const;

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
