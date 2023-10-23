#ifndef Q9_H
#define Q9_H

#include <iostream>
#include <numeric>

#include "Element.hpp"
#include "Node.hpp"
#include "Quadrature.hpp"

// First shear deformation theory second-order quadrangle lagrange element
class Q9 : public Element
{
public:
    // Constructors and destructors //
    Q9(const std::size_t tag) : Element(type_,
                                        tag,
                                        n_nodes_,
                                        dof_per_node_){};

    // The rule of five
    Q9() = default;
    Q9(Q9 const &) = default;
    Q9 &operator=(Q9 const &) = default;
    Q9(Q9 &&) = default;
    Q9 &operator=(Q9 &&) = default;

    // Virtual override //

    // Jacobian transformation methods

    std::vector<double> jacobian_matrix(const double xi,
                                        const double eta) const override;

    double jacobian(const double xi, const double eta) const override;

    std::vector<double> jacobian_inverse_matrix(const double xi,
                                                const double eta) const override;

    // Matrices

    Eigen::MatrixXd stiffness_matrix() const override;
    // Eigen::MatrixXd mass_matrix() const override;

private:
    std::vector<double> shape_functions(const double xi, const double eta) const;
    std::vector<double> shape_functions_d_xi(const double xi,
                                             const double eta) const;
    std::vector<double> shape_functions_d_eta(const double xi,
                                              const double eta) const;
    std::vector<double> shape_functions_d_x(const double xi,
                                            const double eta) const;
    std::vector<double> shape_functions_d_y(const double xi,
                                            const double eta) const;

    // Strain-displacement matrices

    // TODO
    // 1) Test fixed-size Eigen matrices
    // 2) Test Blas/Lapack for storing matrices and performing matrix
    //    multiplication

    // B_m
    Eigen::MatrixXd strain_displacement_membrane(const double xi,
                                                 const double eta) const;

    // B_b
    Eigen::MatrixXd strain_displacement_bending(const double xi,
                                                const double eta) const;

    // B_s
    Eigen::MatrixXd strain_displacement_shear(const double xi,
                                              const double eta) const;

private:
    inline static std::string type_{"Q9"};
    inline static constexpr std::size_t n_nodes_ = 9;
    inline static constexpr std::size_t dof_per_node_ = 5;
};

#endif // Q9_H