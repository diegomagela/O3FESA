#include "Q9.hpp"

// Element shape functions in natural coordinates
std::vector<double> Q9::shape_functions(const double xi, const double eta) const
{

    std::vector<double> shape_function(n_nodes_);

    shape_function.at(0) = (1.0 / 4.0) *
                           eta * xi * (eta * xi - eta - xi + 1.0);

    shape_function.at(1) = (1.0 / 4.0) *
                           eta * xi * (eta * xi + eta - xi - 1.0);

    shape_function.at(2) = (1.0 / 4.0) *
                           eta * xi * (eta * xi + eta + xi + 1.0);

    shape_function.at(3) = (1.0 / 4.0) *
                           eta * xi * (eta * xi - eta + xi - 1.0);

    shape_function.at(4) = (1.0 / 2.0) *
                           eta * (-eta * xi * xi + eta + xi * xi - 1.0);

    shape_function.at(5) = (1.0 / 2.0) *
                           xi * (-eta * eta * xi - eta * eta + xi + 1.0);

    shape_function.at(6) = (1.0 / 2.0) *
                           eta * (-eta * xi * xi + eta - xi * xi + 1.0);

    shape_function.at(7) = (1.0 / 2.0) *
                           xi * (-eta * eta * xi + eta * eta + xi - 1.0);

    shape_function.at(8) = eta * eta * xi * xi - eta * eta - xi * xi + 1.0;

    return shape_function;
}

// Element shape functions derivatives in relation to xi in natural coordinates
std::vector<double> Q9::shape_functions_d_xi(const double xi,
                                             const double eta) const
{
    std::vector<double> shape_function_d_xi(n_nodes_);

    shape_function_d_xi.at(0) = (1.0 / 4.0) *
                                eta * (2.0 * eta * xi - eta - 2.0 * xi + 1.0);

    shape_function_d_xi.at(1) = (1.0 / 4.0) *
                                eta * (2.0 * eta * xi + eta - 2.0 * xi - 1.0);

    shape_function_d_xi.at(2) = (1.0 / 4.0) *
                                eta * (2.0 * eta * xi + eta + 2.0 * xi + 1.0);

    shape_function_d_xi.at(3) = (1.0 / 4.0) *
                                eta * (2.0 * eta * xi - eta + 2.0 * xi - 1.0);

    shape_function_d_xi.at(4) = eta * xi * (1.0 - eta);

    shape_function_d_xi.at(5) = -eta * eta * xi - 1.0 / 2.0 * eta * eta + xi +
                                1.0 / 2.0;

    shape_function_d_xi.at(6) = eta * xi * (-eta - 1.0);

    shape_function_d_xi.at(7) = -eta * eta * xi + (1.0 / 2.0) * eta * eta + xi -
                                1.0 / 2.0;

    shape_function_d_xi.at(8) = 2.0 * xi * (eta * eta - 1.0);

    return shape_function_d_xi;
}

// Element shape functions derivatives in relation to eta in natural coordinates
std::vector<double> Q9::shape_functions_d_eta(const double xi,
                                              const double eta) const
{
    std::vector<double> shape_function_d_eta(n_nodes_);

    shape_function_d_eta.at(0) = (1.0 / 4.0) *
                                 xi * (2.0 * eta * xi - 2.0 * eta - xi + 1.0);

    shape_function_d_eta.at(1) = (1.0 / 4.0) *
                                 xi * (2.0 * eta * xi + 2.0 * eta - xi - 1.0);

    shape_function_d_eta.at(2) = (1.0 / 4.0) *
                                 xi * (2.0 * eta * xi + 2.0 * eta + xi + 1.0);

    shape_function_d_eta.at(3) = (1.0 / 4.0) *
                                 xi * (2.0 * eta * xi - 2.0 * eta + xi - 1.0);

    shape_function_d_eta.at(4) = -eta * xi * xi + eta + (1.0 / 2.0) * xi * xi -
                                 1.0 / 2.0;

    shape_function_d_eta.at(5) = eta * xi * (-xi - 1.0);

    shape_function_d_eta.at(6) = -eta * xi * xi + eta - 1.0 / 2.0 * xi * xi +
                                 1.0 / 2.0;

    shape_function_d_eta.at(7) = eta * xi * (1.0 - xi);

    shape_function_d_eta.at(8) = 2.0 * eta * (xi * xi - 1.0);

    return shape_function_d_eta;
}