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

    std::vector<double>
    jacobian_inverse_matrix(const double xi,
                            const double eta) const override;

    // Matrices
    Eigen::MatrixXd shape_functions_matrix(const double xi,
                                           const double eta) const override;

    Eigen::MatrixXd linear_stiffness_matrix() const override;

    // w_e: the w displacement within the element
    Eigen::MatrixXd
    nonlinear_stiffness_matrix(const std::vector<double> &w_e) const override;

    Eigen::MatrixXd thermal_stiffness_matrix() const;

    // Element tangent matrix
    // q_e: all dofs displacement
    Eigen::MatrixXd
    tangent_stiffness_matrix(const std::vector<double> &q_e) const override;

    // TESTING

    // w_e: the w displacement within the element
    Eigen::VectorXd
    internal_force(const std::vector<double> &w_e) const override;

    // Eigen::MatrixXd mass_matrix() const override;

private:
    std::vector<double> shape_functions(const double xi,
                                        const double eta) const;

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

    // B_theta
    Eigen::MatrixXd strain_displacement_rotation(const double xi,
                                                 const double eta) const;

    // Theta 3x2 matrix
    Eigen::MatrixXd
    strain_displacement_nonlinear_rotation(const double xi,
                                           const double eta,
                                           const std::vector<double> &q_e) const;

    // N_th
    Eigen::VectorXd thermal_force_resultants() const;

    // M_th
    Eigen::VectorXd thermal_moment_resultants() const;

    // F_pressure
    inline Eigen::VectorXd pressure_load_vector() const override;

    // F_thermal
    inline Eigen::VectorXd thermal_load_vector() const override;

private:
    inline static std::string type_{"Q9"};
    inline static constexpr std::size_t n_nodes_ = 9;
    inline static constexpr std::size_t dof_per_node_ = 5;
};

#endif // Q9_H