#include "Q9.hpp"

// Public member functions //

// Private member functions //

// Element shape functions in natural coordinates
std::vector<double> Q9::shape_functions(const double xi, const double eta) const
{
    double N1 = (1.0 / 4.0) * eta * xi * (eta * xi - eta - xi + 1.0);
    double N2 = (1.0 / 4.0) * eta * xi * (eta * xi + eta - xi - 1.0);
    double N3 = (1.0 / 4.0) * eta * xi * (eta * xi + eta + xi + 1.0);
    double N4 = (1.0 / 4.0) * eta * xi * (eta * xi - eta + xi - 1.0);
    double N5 = (1.0 / 2.0) * eta * (-eta * xi * xi + eta + xi * xi - 1.0);
    double N6 = (1.0 / 2.0) * xi * (-eta * eta * xi - eta * eta + xi + 1.0);
    double N7 = (1.0 / 2.0) * eta * (-eta * xi * xi + eta - xi * xi + 1.0);
    double N8 = (1.0 / 2.0) * xi * (-eta * eta * xi + eta * eta + xi - 1.0);
    double N9 = eta * eta * xi * xi - eta * eta - xi * xi + 1.0;

    return {N1, N2, N3, N4, N5, N6, N7, N8, N9};
}

// Element shape functions derivatives in relation to xi in natural coordinates
std::vector<double> Q9::shape_functions_d_xi(const double xi,
                                             const double eta) const
{
    std::vector<double> shape_function_d_xi(n_nodes_);

    double dN1_dxi = (1.0 / 4.0) * eta * (2.0 * eta * xi - eta - 2.0 * xi + 1.0);
    double dN2_dxi = (1.0 / 4.0) * eta * (2.0 * eta * xi + eta - 2.0 * xi - 1.0);
    double dN3_dxi = (1.0 / 4.0) * eta * (2.0 * eta * xi + eta + 2.0 * xi + 1.0);
    double dN4_dxi = (1.0 / 4.0) * eta * (2.0 * eta * xi - eta + 2.0 * xi - 1.0);
    double dN5_dxi = eta * xi * (1.0 - eta);
    double dN6_dxi = -eta * eta * xi - 1.0 / 2.0 * eta * eta + xi + 1.0 / 2.0;
    double dN7_dxi = eta * xi * (-eta - 1.0);
    double dN8_dxi = -eta * eta * xi + (1.0 / 2.0) * eta * eta + xi - 1.0 / 2.0;
    double dN9_dxi = 2.0 * xi * (eta * eta - 1.0);

    return {dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi, dN5_dxi, dN6_dxi, dN7_dxi,
            dN8_dxi, dN9_dxi};
}

// Element shape functions derivatives in relation to eta in natural coordinates
std::vector<double> Q9::shape_functions_d_eta(const double xi,
                                              const double eta) const
{
    double dN1_deta = (1.0 / 4.0) * xi * (2.0 * eta * xi - 2.0 * eta - xi + 1.0);
    double dN2_deta = (1.0 / 4.0) * xi * (2.0 * eta * xi + 2.0 * eta - xi - 1.0);
    double dN3_deta = (1.0 / 4.0) * xi * (2.0 * eta * xi + 2.0 * eta + xi + 1.0);
    double dN4_deta = (1.0 / 4.0) * xi * (2.0 * eta * xi - 2.0 * eta + xi - 1.0);
    double dN5_deta = -eta * xi * xi + eta + (1.0 / 2.0) * xi * xi - 1.0 / 2.0;
    double dN6_deta = eta * xi * (-xi - 1.0);
    double dN7_deta = -eta * xi * xi + eta - 1.0 / 2.0 * xi * xi + 1.0 / 2.0;
    double dN8_deta = eta * xi * (1.0 - xi);
    double dN9_deta = 2.0 * eta * (xi * xi - 1.0);

    return {dN1_deta, dN2_deta, dN3_deta, dN4_deta, dN5_deta, dN6_deta, dN7_deta,
            dN8_deta, dN9_deta};
}

std::vector<double> Q9::shape_functions_d_x(const double xi,
                                            double const eta) const
{
    // |dN_x| = | J_11_inv J_12_inv | |dN_xi |
    // |dN_y|   | J_21_inv J_22_inv | |dN_eta|

    std::vector<double> dN_dxi = shape_functions_d_xi(xi, eta);
    std::vector<double> dN_deta = shape_functions_d_eta(xi, eta);

    double J11_inv = jacobian_inverse_matrix(xi, eta).at(0);
    double J12_inv = jacobian_inverse_matrix(xi, eta).at(1);

    std::vector<double> dN_dx(n_nodes_);

    for (std::size_t i = 0; i < n_nodes_; i++)
        dN_dx.at(i) = J11_inv * dN_dxi.at(i) + J12_inv * dN_deta.at(i);

    return dN_dx;
}

std::vector<double> Q9::shape_functions_d_y(const double xi,
                                            double const eta) const
{
    // |dN_x| = | J_11_inv J_12_inv | |dN_xi |
    // |dN_y|   | J_21_inv J_22_inv | |dN_eta|

    std::vector<double> dN_dxi = shape_functions_d_xi(xi, eta);
    std::vector<double> dN_deta = shape_functions_d_eta(xi, eta);

    double J21_inv = jacobian_inverse_matrix(xi, eta).at(2);
    double J22_inv = jacobian_inverse_matrix(xi, eta).at(3);

    std::vector<double> dN_dy(n_nodes_);

    for (std::size_t i = 0; i < n_nodes_; i++)
        dN_dy.at(i) = J21_inv * dN_dxi.at(i) + J22_inv * dN_deta.at(i);

    return dN_dy;
}

std::vector<double> Q9::jacobian_matrix(const double xi, const double eta) const
{
    std::vector<double> shape_d_xi = shape_functions_d_xi(xi, eta);
    std::vector<double> shape_d_eta = shape_functions_d_eta(xi, eta);

    std::vector<double> x = get_x_coordinates();
    std::vector<double> y = get_y_coordinates();

    double J11 = std::inner_product(x.begin(), x.end(), shape_d_xi.begin(), 0.0);
    double J12 = std::inner_product(y.begin(), y.end(), shape_d_xi.begin(), 0.0);
    double J21 = std::inner_product(x.begin(), x.end(), shape_d_eta.begin(), 0.0);
    double J22 = std::inner_product(y.begin(), y.end(), shape_d_eta.begin(), 0.0);

    return {J11, J12, J21, J22};
}

double Q9::jacobian(const double xi, const double eta) const
{
    double J11 = jacobian_matrix(xi, eta).at(0);
    double J12 = jacobian_matrix(xi, eta).at(1);
    double J21 = jacobian_matrix(xi, eta).at(2);
    double J22 = jacobian_matrix(xi, eta).at(3);

    double J = J11 * J22 - J12 * J21;

    return J;
}

std::vector<double> Q9::jacobian_inverse_matrix(const double xi,
                                                const double eta) const
{
    double J = jacobian(xi, eta);

    double J11 = jacobian_matrix(xi, eta).at(0);
    double J12 = jacobian_matrix(xi, eta).at(1);
    double J21 = jacobian_matrix(xi, eta).at(2);
    double J22 = jacobian_matrix(xi, eta).at(3);

    double J11_inv = (1.0 / J) * J22;
    double J12_inv = (-1.0 / J) * J21;
    double J21_inv = (-1.0 / J) * J12;
    double J22_inv = (1.0 / J) * J11;

    return {J11_inv, J12_inv, J21_inv, J22_inv};
}

Eigen::MatrixXd Q9::strain_displacement_membrane(const double xi,
                                                 const double eta) const
{
    Eigen::VectorXd dN_dx, dN_dy;

    dN_dx = Eigen::Map<Eigen::VectorXd>(
        shape_functions_d_x(xi, eta).data(), n_nodes_);

    dN_dy = Eigen::Map<Eigen::VectorXd>(
        shape_functions_d_y(xi, eta).data(), n_nodes_);

    Eigen::MatrixXd membrane_matrix =
        Eigen::MatrixXd::Zero(3, dof_per_node_ * n_nodes_);

    membrane_matrix.block(0, 0 * n_nodes_, 1, n_nodes_) = dN_dx.transpose();
    membrane_matrix.block(1, 1 * n_nodes_, 1, n_nodes_) = dN_dy.transpose();
    membrane_matrix.block(2, 0 * n_nodes_, 1, n_nodes_) = dN_dy.transpose();
    membrane_matrix.block(2, 1 * n_nodes_, 1, n_nodes_) = dN_dx.transpose();

    return membrane_matrix;
}

Eigen::MatrixXd Q9::strain_displacement_bending(const double xi,
                                                const double eta) const
{
    Eigen::VectorXd dN_dx, dN_dy;

    dN_dx = Eigen::Map<Eigen::VectorXd>(
        shape_functions_d_x(xi, eta).data(), n_nodes_);

    dN_dy = Eigen::Map<Eigen::VectorXd>(
        shape_functions_d_y(xi, eta).data(), n_nodes_);

    Eigen::MatrixXd bending_matrix =
        Eigen::MatrixXd::Zero(3, dof_per_node_ * n_nodes_);

    bending_matrix.block(0, 3 * n_nodes_, 1, n_nodes_) = dN_dx.transpose();
    bending_matrix.block(1, 4 * n_nodes_, 1, n_nodes_) = dN_dy.transpose();
    bending_matrix.block(2, 3 * n_nodes_, 1, n_nodes_) = dN_dy.transpose();
    bending_matrix.block(2, 4 * n_nodes_, 1, n_nodes_) = dN_dx.transpose();

    return bending_matrix;
}

Eigen::MatrixXd Q9::strain_displacement_shear(const double xi,
                                              const double eta) const
{
    Eigen::VectorXd N, dN_dx, dN_dy;

    N = Eigen::Map<Eigen::VectorXd>(
        shape_functions(xi, eta).data(), n_nodes_);

    dN_dx = Eigen::Map<Eigen::VectorXd>(
        shape_functions_d_x(xi, eta).data(), n_nodes_);

    dN_dy = Eigen::Map<Eigen::VectorXd>(
        shape_functions_d_y(xi, eta).data(), n_nodes_);

    Eigen::MatrixXd shear_matrix =
        Eigen::MatrixXd::Zero(2, dof_per_node_ * n_nodes_);

    shear_matrix.block(0, 2 * n_nodes_, 1, n_nodes_) = dN_dx.transpose();
    shear_matrix.block(0, 3 * n_nodes_, 1, n_nodes_) = N.transpose();
    shear_matrix.block(1, 2 * n_nodes_, 1, n_nodes_) = dN_dy.transpose();
    shear_matrix.block(1, 4 * n_nodes_, 1, n_nodes_) = N.transpose();

    return shear_matrix;
}

Eigen::MatrixXd Q9::stiffness_matrix() const
{
    // Material stiffness matrices

    std::vector<double> A_vec = get_section().get()->extensional_stiffness();
    const double A11 = A_vec.at(0);
    const double A12 = A_vec.at(1);
    const double A16 = A_vec.at(2);
    const double A22 = A_vec.at(3);
    const double A26 = A_vec.at(4);
    const double A66 = A_vec.at(5);

    std::vector<double> B_vec = get_section().get()->bending_extensional_stiffness();
    const double B11 = B_vec.at(0);
    const double B12 = B_vec.at(1);
    const double B16 = B_vec.at(2);
    const double B22 = B_vec.at(3);
    const double B26 = B_vec.at(4);
    const double B66 = B_vec.at(5);

    std::vector<double> D_vec = get_section().get()->bending_stiffness();
    const double D11 = D_vec.at(0);
    const double D12 = D_vec.at(1);
    const double D16 = D_vec.at(2);
    const double D22 = D_vec.at(3);
    const double D26 = D_vec.at(4);
    const double D66 = D_vec.at(5);

    std::vector<double> As_vec = get_section().get()->extensional_shear_stiffness();
    // Shear correction factor
    constexpr double Ks = 5.0 / 6.0;

    const double A44 = Ks * As_vec.at(0);
    const double A45 = Ks * As_vec.at(1);
    const double A55 = Ks * As_vec.at(2);

    Eigen::Matrix3d A, B, D;
    Eigen::Matrix2d As;

    A << A11, A12, A16,
        A12, A22, A26,
        A16, A26, A66;

    B << B11, B12, B16,
        B12, B22, B26,
        B16, B26, B66;

    D << D11, D12, D16,
        D12, D22, D26,
        D16, D26, D66;

    As << A44, A45,
        A45, A55;

    constexpr std::size_t n_points = 3;

    Quadrature quadrature(n_points);

    std::vector<double> points = quadrature.get_points();
    std::vector<double> weights = quadrature.get_weights();

    Eigen::MatrixXd K_mm =
        Eigen::MatrixXd::Zero(n_nodes_ * dof_per_node_,
                              n_nodes_ * dof_per_node_);
    Eigen::MatrixXd K_bm =
        Eigen::MatrixXd::Zero(n_nodes_ * dof_per_node_,
                              n_nodes_ * dof_per_node_);
    Eigen::MatrixXd K_bb =
        Eigen::MatrixXd::Zero(n_nodes_ * dof_per_node_,
                              n_nodes_ * dof_per_node_);
    Eigen::MatrixXd K_ss =
        Eigen::MatrixXd::Zero(n_nodes_ * dof_per_node_,
                              n_nodes_ * dof_per_node_);

    for (std::size_t i = 0; i < n_points; ++i)
    {
        for (std::size_t j = 0; j < n_points; ++j)
        {
            double xi_i = points.at(i);
            double eta_j = points.at(j);
            double w_i = weights.at(i);
            double w_j = weights.at(j);

            Eigen::MatrixXd B_m = strain_displacement_membrane(xi_i, eta_j);
            Eigen::MatrixXd B_b = strain_displacement_bending(xi_i, eta_j);
            Eigen::MatrixXd B_s = strain_displacement_shear(xi_i, eta_j);

            // det(J)
            double J = jacobian(xi_i, eta_j);

            // area
            double d_omega = J * w_i * w_j;

            K_mm += B_m.transpose() * A * B_m * d_omega;
            K_bm += B_b.transpose() * B * B_m * d_omega;
            K_bb += B_b.transpose() * D * B_b * d_omega;
            K_ss += B_s.transpose() * As * B_s * d_omega;
        }
    }

    Eigen::MatrixXd K = K_mm + K_bm + K_bm.transpose() + K_bb + K_ss;

    return K;
}