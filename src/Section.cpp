#include "Section.hpp"

#include <iostream>
#include <cmath>
#include <numeric>

// Public functions

inline double Section::total_thickness() const
{
    return std::accumulate(thickness_.begin(), thickness_.end(), 0.0);
}

std::vector<double> Section::layer_position() const
{
    std::vector<double> position(number_of_layers() + 1);

    position.at(0) = -1.0 * total_thickness() / 2.0;

    for (size_t i = 1; i < position.size(); ++i)
        position.at(i) = position.at(i - 1) + total_thickness();

    return position;
}

// Private functions

std::vector<double> Section::transformed_reduced_stiffness(double orientation,
                                                           std::vector<double> elastic_coefficients) const
{
    constexpr double pi = M_PI;

    orientation *= pi / 180.0;

    double Q_11 = elastic_coefficients.at(0);
    double Q_12 = elastic_coefficients.at(1);
    double Q_22 = elastic_coefficients.at(2);
    double Q_66 = elastic_coefficients.at(3);
    double Q_44 = elastic_coefficients.at(4);
    double Q_55 = elastic_coefficients.at(5);

    // c: cos
    // s: sine

    double c_4 = std::pow(std::cos(orientation), 4);
    double c_3 = std::pow(std::cos(orientation), 3);
    double c_2 = std::pow(std::cos(orientation), 2);
    double c = std::cos(orientation);

    double s_4 = std::pow(std::sin(orientation), 4);
    double s_3 = std::pow(std::sin(orientation), 3);
    double s_2 = std::pow(std::sin(orientation), 2);
    double s = std::sin(orientation);

    double Q_bar_11 = Q_11 * c_4 +
                      2.0 * (Q_12 + 2.0 * Q_66) * s_2 * c_2 +
                      Q_22 * s_4;

    double Q_bar_12 = (Q_11 + Q_22 - 4.0 * Q_66) * s_2 * c_2 +
                      Q_12 * (s_4 + c_4);

    double Q_bar_22 = Q_11 * s_4 +
                      2.0 * (Q_12 + 2.0 * Q_66) * s_2 * c_2 +
                      Q_22 * c_4;

    double Q_bar_16 = (Q_11 - Q_12 - 2.0 * Q_66) * s * c_3 +
                      (Q_12 - Q_22 + 2.0 * Q_66) * s_3 * c;

    double Q_bar_26 = (Q_11 - Q_12 - 2.0 * Q_66) * s_3 * c +
                      (Q_12 - Q_22 + 2.0 * Q_66) * s * c_3;

    double Q_bar_66 = (Q_11 + Q_22 - 2.0 * Q_12 - 2.0 * Q_66) * s_2 * c_2 +
                      Q_66 * (s_4 + c_4);

    double Q_bar_44 = Q_44 * c_2 + Q_55 * s_2;
    double Q_bar_45 = (Q_55 - Q_44) * c * s;
    double Q_bar_55 = Q_55 * c_2 + Q_44 * s_2;

    std::vector<double> Q_bar{Q_bar_11, Q_bar_12, Q_bar_22,
                              Q_bar_16, Q_bar_26, Q_bar_66,
                              Q_bar_44, Q_bar_45, Q_bar_55};

    return Q_bar;
}

// Return (symmetric) matrix A
// [A_11, A_12, A_22, A_16, A_26, A_66]
// \sum_{k=1}^N \bar{Q}_{ij} (z_{k+1} - z_{k})
std::vector<double> Section::extensional_stiffness() const
{
    std::size_t size{6};
    std::vector<double> extensional_stiffness(size, 0.0);

    for (std::size_t i = 0; i < number_of_layers(); ++i)
    {
        double orientation = orientation_.at(i);
        double thickness = thickness_.at(i);
        auto elastic_coefficients = material_.at(i).get()->elastic_coefficients();

        for (std::size_t i = 0; i < size; i++)
        {
            extensional_stiffness.at(i) +=
                transformed_reduced_stiffness(orientation, elastic_coefficients).at(i) * thickness;
        }
    }

    return extensional_stiffness;
}

// Return (symmetric) matrix D
// [D_11, D_12, D_22, D_16, D_26, D_66]
// \sum_{k=1}^N \bar{Q}_{ij} (z^3_{k+1} - z^3_{k})
std::vector<double> Section::bending_stiffness() const
{
    std::size_t size{6};
    std::vector<double> bending_stiffness(size, 0.0);

    for (std::size_t i = 0; i < number_of_layers(); ++i)
    {
        double orientation = orientation_.at(i);
        double z_1 = layer_position().at(i + 1);
        double z_0 = layer_position().at(i + 1);
        double z_diff = std::pow(z_1, 3) - std::pow(z_0, 3);

        auto elastic_coefficients = material_.at(i).get()->elastic_coefficients();

        for (std::size_t i = 0; i < size; i++)
        {
            bending_stiffness.at(i) += (1.0 / 3.0) *
                                       transformed_reduced_stiffness(orientation, elastic_coefficients).at(i) * z_diff;
        }
    }

    return bending_stiffness;
}

// Return (symmetric) matrix B
// [B_11, B_12, B_22, B_16, B_26, B_66]
// \sum_{k=1}^N \bar{Q}_{ij} (z^2_{k+1} - z^2_{k})
std::vector<double> Section::bending_extensional_stiffness() const
{
    std::size_t size{6};
    std::vector<double> bending_extensional_stiffness(size, 0.0);

    for (std::size_t i = 0; i < number_of_layers(); ++i)
    {
        double orientation = orientation_.at(i);
        double z_1 = layer_position().at(i + 1);
        double z_0 = layer_position().at(i + 1);
        double z_diff = std::pow(z_1, 2) - std::pow(z_0, 2);

        auto elastic_coefficients = material_.at(i).get()->elastic_coefficients();

        for (std::size_t i = 0; i < size; i++)
        {
            bending_extensional_stiffness.at(i) += (1.0 / 2.0) *
                                                   transformed_reduced_stiffness(orientation, elastic_coefficients).at(i) * z_diff;
        }
    }

    return bending_extensional_stiffness;
}

// Return (symmetric) matrix C
// [A_44, A_45, A_55]
// \sum_{k=1}^N \bar{Q}_{ij} (z_{k+1} - z_{k})
std::vector<double> Section::extensional_shear_stiffness() const
{
    std::size_t size{3};
    std::vector<double> extensional_shear_stiffness(size, 0.0);

    for (std::size_t i = 0; i < number_of_layers(); ++i)
    {
        double orientation = orientation_.at(i);
        double thickness = thickness_.at(i);
        auto elastic_coefficients = material_.at(i).get()->elastic_coefficients();

        for (std::size_t i = 0; i < size; i++)
        {
            extensional_shear_stiffness.at(i) +=
                transformed_reduced_stiffness(orientation, elastic_coefficients).at(i) * thickness;
        }
    }

    return extensional_shear_stiffness;
}