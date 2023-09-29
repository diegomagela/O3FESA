#ifndef NODE_H
#define NODE_H

#include "Boundary.hpp"
#include "CLoad.hpp"

#include <cstddef>
#include <iostream>

// Node class owns a node tag and its 3D coordinates
class Node
{
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
    inline bool has_boundary() const { return boundary_.has_boundary(); }
    inline bool has_cload() const { return cload_.has_cload(); }
    void print() const;

    // Modifiers
    inline void set_boundary(Boundary boundary) { boundary_ = boundary; };
    inline void set_cload(CLoad cload) { cload_ = cload; }

    // Friends
    friend std::ostream &operator<<(std::ostream &os, const Node &node);

private:
    const std::size_t tag_{};
    double x_{};
    double y_{};
    double z_{};
    Boundary boundary_;
    CLoad cload_;
};

#endif // NODE_H
