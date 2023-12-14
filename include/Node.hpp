#ifndef NODE_H
#define NODE_H

#include "Boundary.hpp"
#include "CLoad.hpp"
#include "Dofs.hpp"

#include <cstddef>
#include <iostream>
#include <memory>

// Node class owns a tag, its 3D coordinates, boundary condition and
// concentrated loading
class Node
{
private:
    typedef std::shared_ptr<Boundary> BoundaryPtr;
    typedef std::shared_ptr<CLoad> CLoadPtr;

public:
    explicit Node(const std::size_t tag,
                  const double x,
                  const double y,
                  const double z) noexcept : tag_(tag),
                                             x_(x),
                                             y_(y),
                                             z_(z){};

    // The rule of five

    Node() = default;
    Node(Node const &) = default;
    Node &operator=(Node const &) = default;
    Node(Node &&) = default;
    Node &operator=(Node &&) = default;

    // Selectors

    // Returns node's tag
    inline std::size_t node_tag() const { return tag_; };
    // Returns number of degrees of freedom
    inline std::size_t n_dof() const { return n_dof_; }
    // All dofs displacement at each node
    inline std::vector<double>
    get_displacements() const { return dofs_.get_displacements(); }
    // Node's u displacement
    inline double 
    get_u_displacements() const { return dofs_.get_u_displacement(); }
    // Node's v displacement
    inline double 
    get_v_displacements() const { return dofs_.get_v_displacement(); }
    // Node's w displacement
    inline double 
    get_w_displacements() const { return dofs_.get_w_displacement(); }
    // Node's phi_x displacement
    inline double 
    get_phi_x_displacements() const { return dofs_.get_phi_x_displacement(); }
    // Node's phi_y displacement
    inline double 
    get_phi_y_displacements() const { return dofs_.get_phi_y_displacement(); }
    // Returns x coordinate
    inline double get_x() const { return x_; }
    // Returns y coordinate
    inline double get_y() const { return y_; }
    // Returns z coordinate
    inline double get_z() const { return z_; }
    // Returns Boundary class pointer
    inline BoundaryPtr get_boundary() const { return boundary_; }
    // Returns CLoad class pointer
    inline CLoadPtr get_cload() const { return cload_; }
    // Returns if node has boundary dof
    bool has_boundary() const;
    // Returns if node has concentrated load
    bool has_cload() const;
    // Print node data
    void print() const;

    // Modifiers

    // Set Boundary class pointer to node
    inline void set_boundary(BoundaryPtr boundary) { boundary_ = boundary; };
    // Set CLoad class pointer to node
    inline void set_cload(CLoadPtr cload) { cload_ = cload; }
    // Set Dofs class to node
    inline void set_dofs(Dofs dofs) { dofs_ = dofs; }

    // Friends

    // Friend operator<< for printing data
    friend std::ostream &operator<<(std::ostream &os, const Node &node);

private:
    // TODO
    // How to identify how many dofs the node has?
    // Considering 5 DOF per node
    const std::size_t n_dof_ = 5;

    const std::size_t tag_{};
    double x_{};
    double y_{};
    double z_{};
    BoundaryPtr boundary_;
    CLoadPtr cload_;
    Dofs dofs_;
};

#endif // NODE_H
