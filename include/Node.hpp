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
    inline std::size_t node_tag() const { return tag_; };
    inline std::size_t n_dof() const { return n_dof_; }
    
    inline std::vector<double>
    get_displacements() const { return dofs_.get_displacements(); }

    inline double
    get_w_displacements() const { return dofs_.get_w_displacement(); }

    // Return node coordinates
    inline double get_x() const { return x_; }
    inline double get_y() const { return y_; }
    inline double get_z() const { return z_; }
    inline BoundaryPtr get_boundary() const { return boundary_; }
    inline CLoadPtr get_cload() const { return cload_; }

    bool has_boundary() const;
    bool has_cload() const;
    void print() const;

    // Modifiers
    inline void set_boundary(BoundaryPtr boundary) { boundary_ = boundary; };
    inline void set_cload(CLoadPtr cload) { cload_ = cload; }
    inline void set_dofs(Dofs dofs) { dofs_ = dofs; }

    // Friends
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
