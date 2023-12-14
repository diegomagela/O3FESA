#include "Q9.hpp"

#include <fstream>

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

Eigen::MatrixXd Q9::shape_functions_matrix(const double xi,
                                           const double eta) const
{
    Eigen::VectorXd N = vct::std_to_eigen_vector(shape_functions(xi, eta));

    Eigen::MatrixXd shape_functions_matrix =
        Eigen::MatrixXd::Zero(dof_per_node_, dof_per_node_ * n_nodes_);

    shape_functions_matrix.block(0, 0 * n_nodes_, 1, n_nodes_) = N.transpose();
    shape_functions_matrix.block(1, 1 * n_nodes_, 1, n_nodes_) = N.transpose();
    shape_functions_matrix.block(2, 2 * n_nodes_, 1, n_nodes_) = N.transpose();
    shape_functions_matrix.block(3, 3 * n_nodes_, 1, n_nodes_) = N.transpose();
    shape_functions_matrix.block(4, 4 * n_nodes_, 1, n_nodes_) = N.transpose();

    return shape_functions_matrix;
}

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

    return {dN1_deta, dN2_deta, dN3_deta,
            dN4_deta, dN5_deta, dN6_deta,
            dN7_deta, dN8_deta, dN9_deta};
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

    for (std::size_t i = 0; i < n_nodes_; ++i)
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

    for (std::size_t i = 0; i < n_nodes_; ++i)
        dN_dy.at(i) = J21_inv * dN_dxi.at(i) + J22_inv * dN_deta.at(i);

    return dN_dy;
}

std::vector<double> Q9::jacobian_matrix(const double xi, const double eta) const
{
    std::vector<double> shape_d_xi = shape_functions_d_xi(xi, eta);
    std::vector<double> shape_d_eta = shape_functions_d_eta(xi, eta);

    std::vector<double> x = get_x_coordinates();
    std::vector<double> y = get_y_coordinates();

    double J11 =
        std::inner_product(x.begin(), x.end(), shape_d_xi.begin(), 0.0);

    double J12 =
        std::inner_product(y.begin(), y.end(), shape_d_xi.begin(), 0.0);

    double J21 =
        std::inner_product(x.begin(), x.end(), shape_d_eta.begin(), 0.0);

    double J22 =
        std::inner_product(y.begin(), y.end(), shape_d_eta.begin(), 0.0);

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

    dN_dx = vct::std_to_eigen_vector(shape_functions_d_x(xi, eta));

    dN_dy = vct::std_to_eigen_vector(shape_functions_d_y(xi, eta));

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

    dN_dx = vct::std_to_eigen_vector(shape_functions_d_x(xi, eta));

    dN_dy = vct::std_to_eigen_vector(shape_functions_d_y(xi, eta));

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

    N = vct::std_to_eigen_vector(shape_functions(xi, eta));

    dN_dx = vct::std_to_eigen_vector(shape_functions_d_x(xi, eta));

    dN_dy = vct::std_to_eigen_vector(shape_functions_d_y(xi, eta));

    Eigen::MatrixXd shear_matrix =
        Eigen::MatrixXd::Zero(2, dof_per_node_ * n_nodes_);

    shear_matrix.block(0, 2 * n_nodes_, 1, n_nodes_) = dN_dx.transpose();
    shear_matrix.block(0, 3 * n_nodes_, 1, n_nodes_) = N.transpose();
    shear_matrix.block(1, 2 * n_nodes_, 1, n_nodes_) = dN_dy.transpose();
    shear_matrix.block(1, 4 * n_nodes_, 1, n_nodes_) = N.transpose();

    return shear_matrix;
}

Eigen::MatrixXd Q9::strain_displacement_rotation(const double xi,
                                                 const double eta) const
{
    Eigen::VectorXd dN_dx, dN_dy;

    dN_dx = vct::std_to_eigen_vector(shape_functions_d_x(xi, eta));

    dN_dy = vct::std_to_eigen_vector(shape_functions_d_y(xi, eta));

    Eigen::MatrixXd rotation_matrix =
        Eigen::MatrixXd::Zero(2, dof_per_node_ * n_nodes_);

    rotation_matrix.block(0, 2 * n_nodes_, 1, n_nodes_) = dN_dx.transpose();
    rotation_matrix.block(1, 2 * n_nodes_, 1, n_nodes_) = dN_dy.transpose();

    return rotation_matrix;
}

Eigen::MatrixXd
Q9::strain_displacement_nonlinear_rotation(const double xi,
                                           const double eta) const
{
    // Theta =  |dw_dx      0    |
    //          |0          dw_dy|
    //          |dw_dy      dw_dx|

    // w = [N_1 N_2 ... N_m] [w_1 w_2 ... w_m]^T

    Eigen::VectorXd dN_dx, dN_dy, w_e;

    dN_dx = vct::std_to_eigen_vector(shape_functions_d_x(xi, eta));
    dN_dy = vct::std_to_eigen_vector(shape_functions_d_y(xi, eta));

    std::vector<double> w_element = get_w_displacement();

    w_e = vct::std_to_eigen_vector(w_element);

    double dw_dx = dN_dx.dot(w_e);
    double dw_dy = dN_dy.dot(w_e);

    Eigen::MatrixXd theta_matrix = Eigen::MatrixXd::Zero(3, 2);
    theta_matrix(0, 0) = dw_dx;
    theta_matrix(1, 1) = dw_dy;
    theta_matrix(2, 0) = dw_dy;
    theta_matrix(2, 1) = dw_dx;

    return theta_matrix;
}

Eigen::VectorXd Q9::thermal_force_resultants() const
{
    std::vector<double> in_plane_thermal =
        get_section().get()->thermal_force_resultants();

    Eigen::VectorXd N_th =
        vct::std_to_eigen_vector(in_plane_thermal);

    return N_th;
}

Eigen::VectorXd Q9::thermal_moment_resultants() const
{
    std::vector<double> moment_thermal =
        get_section().get()->thermal_moment_resultants();

    Eigen::VectorXd M_th =
        vct::std_to_eigen_vector(moment_thermal);

    return M_th;
}

Eigen::VectorXd Q9::in_plane_force_resultants(const double xi,
                                              const double eta) const
{
    // Material stiffness matrices

    std::vector<double> A_vec = get_section().get()->extensional_stiffness();
    const double A11 = A_vec.at(0);
    const double A12 = A_vec.at(1);
    const double A16 = A_vec.at(2);
    const double A22 = A_vec.at(3);
    const double A26 = A_vec.at(4);
    const double A66 = A_vec.at(5);

    std::vector<double> B_vec =
        get_section().get()->bending_extensional_stiffness();

    const double B11 = B_vec.at(0);
    const double B12 = B_vec.at(1);
    const double B16 = B_vec.at(2);
    const double B22 = B_vec.at(3);
    const double B26 = B_vec.at(4);
    const double B66 = B_vec.at(5);

    Eigen::Matrix3d A, B;

    A << A11, A12, A16,
        A12, A22, A26,
        A16, A26, A66;

    B << B11, B12, B16,
        B12, B22, B26,
        B16, B26, B66;

    // Element's displacements
    Eigen::VectorXd u_element = vct::std_to_eigen_vector(get_u_displacement());
    Eigen::VectorXd v_element = vct::std_to_eigen_vector(get_v_displacement());
    Eigen::VectorXd w_element = vct::std_to_eigen_vector(get_w_displacement());
    Eigen::VectorXd phi_x_element = vct::std_to_eigen_vector(get_phi_x_displacement());
    Eigen::VectorXd phi_y_element = vct::std_to_eigen_vector(get_phi_y_displacement());

    // Shape functions derivatives w.r.t. cartesian coordinates
    Eigen::VectorXd dN_dx = vct::std_to_eigen_vector(shape_functions_d_x(xi, eta));
    Eigen::VectorXd dN_dy = vct::std_to_eigen_vector(shape_functions_d_y(xi, eta));

    // Displacement's derivatives
    double du_dx = dN_dx.dot(u_element);
    double du_dy = dN_dy.dot(u_element);

    double dv_dx = dN_dx.dot(v_element);
    double dv_dy = dN_dy.dot(v_element);

    double dw_dx = dN_dx.dot(w_element);
    double dw_dy = dN_dy.dot(w_element);

    double dphi_x_dx = dN_dx.dot(phi_x_element);
    double dphi_x_dy = dN_dy.dot(phi_x_element);

    double dphi_y_dx = dN_dx.dot(phi_y_element);
    double dphi_y_dy = dN_dy.dot(phi_y_element);

    // Membrane strain \varepsilon_m
    Eigen::Vector3d membrane_strain;
    membrane_strain(0) = du_dx + (1.0 / 2.0) * dw_dx * dw_dx;
    membrane_strain(1) = dv_dy + (1.0 / 2.0) * dw_dy * dw_dy;
    membrane_strain(2) = (du_dy + dv_dx) + dw_dx * dw_dy;

    // Bending strain \varepsilon_b
    Eigen::Vector3d bending_strain;
    bending_strain(0) = dphi_x_dx;
    bending_strain(1) = dphi_y_dy;
    bending_strain(2) = dphi_x_dy + dphi_y_dx;

    // Thermal in-plane force resultants
    Eigen::VectorXd N_T = thermal_force_resultants();

    // In-plane force resultants
    Eigen::VectorXd in_plane = A * membrane_strain + B * bending_strain - N_T;

    return in_plane;
}

Eigen::MatrixXd Q9::linear_stiffness_matrix() const
{
    // Material stiffness matrices

    std::vector<double> A_vec = get_section().get()->extensional_stiffness();
    const double A11 = A_vec.at(0);
    const double A12 = A_vec.at(1);
    const double A16 = A_vec.at(2);
    const double A22 = A_vec.at(3);
    const double A26 = A_vec.at(4);
    const double A66 = A_vec.at(5);

    std::vector<double> B_vec =
        get_section().get()->bending_extensional_stiffness();

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

    constexpr std::size_t n_quadrature_points{3};

    Quadrature quadrature(n_quadrature_points);

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

    for (std::size_t i = 0; i < n_quadrature_points; ++i)
    {
        for (std::size_t j = 0; j < n_quadrature_points; ++j)
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

Eigen::MatrixXd Q9::nonlinear_stiffness_matrix() const
{
    // Element's w displacement
    std::vector<double> w_e = get_w_displacement();

    // Material stiffness matrices

    std::vector<double> A_vec = get_section().get()->extensional_stiffness();
    const double A11 = A_vec.at(0);
    const double A12 = A_vec.at(1);
    const double A16 = A_vec.at(2);
    const double A22 = A_vec.at(3);
    const double A26 = A_vec.at(4);
    const double A66 = A_vec.at(5);

    std::vector<double> B_vec =
        get_section().get()->bending_extensional_stiffness();

    const double B11 = B_vec.at(0);
    const double B12 = B_vec.at(1);
    const double B16 = B_vec.at(2);
    const double B22 = B_vec.at(3);
    const double B26 = B_vec.at(4);
    const double B66 = B_vec.at(5);

    Eigen::Matrix3d A, B;

    A << A11, A12, A16,
        A12, A22, A26,
        A16, A26, A66;

    B << B11, B12, B16,
        B12, B22, B26,
        B16, B26, B66;

    constexpr std::size_t n_quadrature_points{3};

    Quadrature quadrature(n_quadrature_points);

    std::vector<double> points = quadrature.get_points();
    std::vector<double> weights = quadrature.get_weights();

    Eigen::MatrixXd K_mt =
        Eigen::MatrixXd::Zero(n_nodes_ * dof_per_node_,
                              n_nodes_ * dof_per_node_);
    Eigen::MatrixXd K_tt =
        Eigen::MatrixXd::Zero(n_nodes_ * dof_per_node_,
                              n_nodes_ * dof_per_node_);
    Eigen::MatrixXd K_tb =
        Eigen::MatrixXd::Zero(n_nodes_ * dof_per_node_,
                              n_nodes_ * dof_per_node_);

    for (std::size_t i = 0; i < n_quadrature_points; ++i)
    {
        for (std::size_t j = 0; j < n_quadrature_points; ++j)
        {
            double xi_i = points.at(i);
            double eta_j = points.at(j);
            double w_i = weights.at(i);
            double w_j = weights.at(j);

            Eigen::MatrixXd B_m = strain_displacement_membrane(xi_i, eta_j);
            Eigen::MatrixXd B_b = strain_displacement_bending(xi_i, eta_j);
            Eigen::MatrixXd B_t = strain_displacement_rotation(xi_i, eta_j);

            Eigen::MatrixXd Theta =
                strain_displacement_nonlinear_rotation(xi_i, eta_j);

            // det(J)
            double J = jacobian(xi_i, eta_j);

            // area
            double d_omega = J * w_i * w_j;

            K_mt += B_m.transpose() * A * Theta * B_t * d_omega;
            K_tt += B_t.transpose() * Theta.transpose() * A * Theta * B_t *
                    d_omega;
            K_tb += B_t.transpose() * Theta.transpose() * B * B_b * d_omega;
        }
    }

    Eigen::MatrixXd K_NL = (1.0 / 2.0) * (K_mt + K_tt + K_tb.transpose()) +
                           K_tb + K_mt.transpose();

    return K_NL;
}

Eigen::MatrixXd Q9::thermal_stiffness_matrix() const
{
    double N_x = thermal_force_resultants()[0];
    double N_y = thermal_force_resultants()[1];
    double N_xy = thermal_force_resultants()[2];

    // N_{delta T}
    Eigen::MatrixXd N_T = Eigen::MatrixXd::Zero(2, 2);
    N_T(0, 0) = N_x;
    N_T(0, 1) = N_xy;
    N_T(1, 0) = N_xy;
    N_T(1, 1) = N_y;

    constexpr std::size_t n_quadrature_points{3};

    Quadrature quadrature(n_quadrature_points);

    std::vector<double> points = quadrature.get_points();
    std::vector<double> weights = quadrature.get_weights();

    Eigen::MatrixXd K_dT =
        Eigen::MatrixXd::Zero(n_nodes_ * dof_per_node_,
                              n_nodes_ * dof_per_node_);

    for (std::size_t i = 0; i < n_quadrature_points; ++i)
    {
        for (std::size_t j = 0; j < n_quadrature_points; ++j)
        {
            double xi_i = points.at(i);
            double eta_j = points.at(j);
            double w_i = weights.at(i);
            double w_j = weights.at(j);

            Eigen::MatrixXd B_t = strain_displacement_rotation(xi_i, eta_j);

            // det(J)
            double J = jacobian(xi_i, eta_j);

            // area
            double d_omega = J * w_i * w_j;

            K_dT += B_t.transpose() * N_T * B_t * d_omega;
        }
    }

    return -1.0 * K_dT;
}

Eigen::MatrixXd Q9::tangent_stiffness_matrix() const
{
    // Material stiffness matrices

    std::vector<double> A_vec = get_section().get()->extensional_stiffness();
    const double A11 = A_vec.at(0);
    const double A12 = A_vec.at(1);
    const double A16 = A_vec.at(2);
    const double A22 = A_vec.at(3);
    const double A26 = A_vec.at(4);
    const double A66 = A_vec.at(5);

    std::vector<double> B_vec =
        get_section().get()->bending_extensional_stiffness();
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

    Eigen::Matrix3d A, B, D;

    A << A11, A12, A16,
        A12, A22, A26,
        A16, A26, A66;

    B << B11, B12, B16,
        B12, B22, B26,
        B16, B26, B66;

    D << D11, D12, D16,
        D12, D22, D26,
        D16, D26, D66;

    // All elements' dofs displacement
    std::vector<double> q_e = get_displacements();

    // Elements' w displacement
    std::vector<double> w_e = get_w_displacement();

    constexpr std::size_t n_quadrature_points{3};

    Quadrature quadrature(n_quadrature_points);

    std::vector<double> points = quadrature.get_points();
    std::vector<double> weights = quadrature.get_weights();

    Eigen::MatrixXd K_mt = Eigen::MatrixXd::Zero(n_nodes_ * dof_per_node_,
                                                 n_nodes_ * dof_per_node_);

    Eigen::MatrixXd K_tt = Eigen::MatrixXd::Zero(n_nodes_ * dof_per_node_,
                                                 n_nodes_ * dof_per_node_);

    Eigen::MatrixXd K_tb = Eigen::MatrixXd::Zero(n_nodes_ * dof_per_node_,
                                                 n_nodes_ * dof_per_node_);

    // In-plane resultant vector
    Eigen::VectorXd N_vector = Eigen::VectorXd::Zero(n_nodes_ * dof_per_node_);

    // Initial stress matrix
    // K_sigma
    Eigen::MatrixXd K_sg = Eigen::MatrixXd::Zero(n_nodes_ * dof_per_node_,
                                                 n_nodes_ * dof_per_node_);

    Eigen::VectorXd d_e = vct::std_to_eigen_vector(q_e);

    for (std::size_t i = 0; i < n_quadrature_points; ++i)
    {
        for (std::size_t j = 0; j < n_quadrature_points; ++j)
        {
            double xi_i = points.at(i);
            double eta_j = points.at(j);

            double w_i = weights.at(i);
            double w_j = weights.at(j);

            Eigen::MatrixXd B_m = strain_displacement_membrane(xi_i, eta_j);
            Eigen::MatrixXd B_b = strain_displacement_bending(xi_i, eta_j);
            Eigen::MatrixXd B_t = strain_displacement_rotation(xi_i, eta_j);
            Eigen::MatrixXd Theta =
                strain_displacement_nonlinear_rotation(xi_i, eta_j);
            Eigen::VectorXd N_vector = in_plane_force_resultants(xi_i, eta_j);

            // Stress matrix \sigma
            double N_xx = N_vector(0);
            double N_yy = N_vector(1);
            double N_xy = N_vector(2);

            Eigen::Matrix2d N_sg;
            N_sg << N_xx, N_xy, N_xy, N_yy;

            // det(J)
            double J = jacobian(xi_i, eta_j);

            // area
            double d_omega = J * w_i * w_j;

            K_mt += B_m.transpose() * A * Theta * B_t * d_omega;

            K_tt += B_t.transpose() * Theta.transpose() * A * Theta * B_t *
                    d_omega;

            K_tb += B_t.transpose() * Theta.transpose() * B * B_b * d_omega;

            // >>>
            K_sg += B_t.transpose() * N_sg * B_t * d_omega;
        }
    }

    Eigen::MatrixXd T = K_mt + K_mt.transpose() +
                        K_tb + K_tb.transpose() +
                        K_tt +
                        K_sg;

    return T;
}

Eigen::VectorXd Q9::internal_force_vector() const
{
    Eigen::MatrixXd nonlinear_matrix = nonlinear_stiffness_matrix();
    Eigen::MatrixXd linear_matrix = linear_stiffness_matrix();

    std::vector<double> element_displacement = get_displacements();
    Eigen::VectorXd q_e = vct::std_to_eigen_vector(element_displacement);

    // As the residual is defined as 
    // R = f_int - f_ext
    // and
    // K_T \delta q = -R
    // f_int is already computed as -f_int

    Eigen::VectorXd internal_force = -1.0 * (nonlinear_matrix + linear_matrix) * q_e;

    return internal_force;
}

Eigen::VectorXd Q9::pressure_load_vector() const
{
    Eigen::VectorXd load_vector =
        Eigen::VectorXd::Zero(dof_per_node_ * n_nodes_);

    if (has_load())
    {
        constexpr std::size_t n_quadrature_points = 2;

        Quadrature quadrature(n_quadrature_points);

        std::vector<double> points = quadrature.get_points();
        std::vector<double> weights = quadrature.get_weights();

        // Loading vector
        // f^e = N^T p
        // It is being considered only transversal distributed loading

        const double p_w = load_value();

        Eigen::VectorXd p(dof_per_node_);
        p << 0, 0, p_w, 0, 0;

        for (std::size_t i = 0; i < n_quadrature_points; ++i)
        {
            for (std::size_t j = 0; j < n_quadrature_points; ++j)
            {
                double xi_i = points.at(i);
                double eta_j = points.at(j);
                double w_i = weights.at(i);
                double w_j = weights.at(j);

                Eigen::MatrixXd N = shape_functions_matrix(xi_i, eta_j);

                // det(J)
                double J = jacobian(xi_i, eta_j);

                // area
                double d_omega = J * w_i * w_j;

                load_vector += N.transpose() * p * d_omega;
            }
        }
    }

    return load_vector;
}

Eigen::VectorXd Q9::thermal_load_vector() const
{
    Eigen::VectorXd load_vector =
        Eigen::VectorXd::Zero(dof_per_node_ * n_nodes_);

    if (has_thermal_load())
    {
        constexpr std::size_t n_quadrature_points{2};

        Quadrature quadrature(n_quadrature_points);

        std::vector<double> points = quadrature.get_points();
        std::vector<double> weights = quadrature.get_weights();

        Eigen::VectorXd N_th = thermal_force_resultants();
        Eigen::VectorXd M_th = thermal_moment_resultants();

        for (std::size_t i = 0; i < n_quadrature_points; ++i)
        {
            for (std::size_t j = 0; j < n_quadrature_points; ++j)
            {
                double xi_i = points.at(i);
                double eta_j = points.at(j);
                double w_i = weights.at(i);
                double w_j = weights.at(j);

                Eigen::MatrixXd Dm = strain_displacement_membrane(xi_i, eta_j);
                Eigen::MatrixXd Db = strain_displacement_bending(xi_i, eta_j);

                // det(J)
                double J = jacobian(xi_i, eta_j);

                // area
                double d_omega = J * w_i * w_j;

                load_vector += (Dm.transpose() * N_th + Db.transpose() * M_th) *
                               d_omega;
            }
        }
    }

    return load_vector;
}

Eigen::VectorXd Q9::external_load_vector() const
{
    return pressure_load_vector() + thermal_load_vector();
}

std::vector<Triplet> Q9::linear_stiffness_triplet() const
{
    return matrix_to_triplet(linear_stiffness_matrix());
}

std::vector<Triplet> Q9::nonlinear_stiffness_triplet() const
{
    return matrix_to_triplet(nonlinear_stiffness_matrix());
}

std::vector<Triplet> Q9::external_load_triplet() const
{
    return vector_to_triplet(external_load_vector());
}

std::vector<Triplet> Q9::internal_load_triplet() const
{
    return vector_to_triplet(internal_force_vector());
}

std::vector<Triplet> Q9::tangent_stiffness_triplet() const
{
    return matrix_to_triplet(tangent_stiffness_matrix());
}