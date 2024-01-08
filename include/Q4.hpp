#ifndef Q4_H
#define Q4_H

#include "Element.hpp"
#include "Node.hpp"
#include "Quadrature.hpp"

#include <iostream>
#include <numeric>

/*
    TODO
    - Test fixed-size Eigen matrices
    - Test Blas for performing linear algebra operations
*/

// 4-node lagrange element
class Q4 : public Element
{
public:
    // Constructors and destructors //
    Q4(const std::size_t tag) : Element(type_,
                                        tag,
                                        n_nodes_,
                                        dof_per_node_){};

    // The rule of five
    Q4() = default;
    Q4(Q4 const &) = default;
    Q4 &operator=(Q4 const &) = default;
    Q4(Q4 &&) = default;
    Q4 &operator=(Q4 &&) = default;

private:
    // Q4 shape functions
    std::vector<double> shape_functions(const double xi,
                                        const double eta) const;
    // Matrix containing all dofs shape function
    Eigen::MatrixXd shape_functions_matrix(const double xi,
                                           const double eta) const;
    // Shape functions derivative w.r.t. natural coordinate \\xi
    std::vector<double> shape_functions_d_xi(const double xi,
                                             const double eta) const;
    // Shape functions derivative w.r.t. natural coordinate \\eta
    std::vector<double> shape_functions_d_eta(const double xi,
                                              const double eta) const;
    // Shape functions derivative w.r.t. cartesian coordinate x
    std::vector<double> shape_functions_d_x(const double xi,
                                            const double eta) const;
    // Shape functions derivative w.r.t. cartesian coordinate y
    std::vector<double> shape_functions_d_y(const double xi,
                                            const double eta) const;
    // B_m matrix
    Eigen::MatrixXd strain_displacement_membrane(const double xi,
                                                 const double eta) const;
    // B_b matrix
    Eigen::MatrixXd strain_displacement_bending(const double xi,
                                                const double eta) const;
    // B_s matrix
    Eigen::MatrixXd strain_displacement_shear(const double xi,
                                              const double eta) const;
    // B_\\theta matrix
    Eigen::MatrixXd strain_displacement_rotation(const double xi,
                                                 const double eta) const;
    // \\Theta matrix
    Eigen::MatrixXd
    strain_displacement_nonlinear_rotation(const double xi,
                                           const double eta) const;
    // In-plane resultant forces
    Eigen::VectorXd in_plane_force_resultants(const double xi,
                                              const double eta) const;
    // Jacobian transformation matrix
    std::vector<double> jacobian_matrix(const double xi,
                                        const double eta) const;
    // Jacobian matrix determinant
    double jacobian(const double xi, const double eta) const;
    // Jacobian transverse inverse matrix
    std::vector<double> jacobian_inverse_matrix(const double xi,
                                                const double eta) const;
    // Linear dense stiffness matrix
    Eigen::MatrixXd linear_stiffness_matrix() const;
    // Nonlinear dense stiffness matrix (nonlinear part)
    Eigen::MatrixXd nonlinear_stiffness_matrix() const;
    // Thermal dense stiffness matrix
    Eigen::MatrixXd thermal_stiffness_matrix() const;
    // Internal dense force vector K(q)*q
    Eigen::VectorXd internal_force_vector() const;
    // Tangent dense stiffness matrix (nonlinear part)
    Eigen::MatrixXd tangent_stiffness_matrix() const;
    // Resultant in-plane thermal forces
    Eigen::VectorXd thermal_force_resultants() const;
    // Resultant thermal moments
    Eigen::VectorXd thermal_moment_resultants() const;
    // Distributed pressure load
    Eigen::VectorXd pressure_load_vector() const;
    // Element thermal load
    Eigen::VectorXd thermal_load_vector() const;
    // External load
    Eigen::VectorXd external_load_vector() const;

    // Virtual functions override

    // Triplet linear matrix
    std::vector<Triplet> linear_stiffness_triplet() const override;
    // Triplet nonlinear matrix
    std::vector<Triplet> nonlinear_stiffness_triplet() const override;
    // Triplet external load vector
    std::vector<Triplet> external_load_triplet() const override;
    // Triplet internal load vector: K(q)q
    std::vector<Triplet> internal_load_triplet() const override;
    // Triplet tangent stiffness matrix
    std::vector<Triplet> tangent_stiffness_triplet() const override;

private:
    inline static std::string type_{"Q4"};
    inline static constexpr std::size_t n_nodes_ = 4;
    inline static constexpr std::size_t dof_per_node_ = 5;
};

#endif // Q9_H