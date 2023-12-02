#ifndef NODE_H
#define NODE_H

#include "Boundary.hpp"
#include "CLoad.hpp"

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

    // Friends
    friend std::ostream &operator<<(std::ostream &os, const Node &node);

private:
    const std::size_t tag_{};
    double x_{};
    double y_{};
    double z_{};
    BoundaryPtr boundary_;
    CLoadPtr cload_;
};

#endif // NODE_H
