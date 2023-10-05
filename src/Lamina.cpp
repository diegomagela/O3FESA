#include "Lamina.hpp"

std::vector<double> Lamina::elastic_coefficients() const
{
    // E_1, E_2, NU_12, G_12, G_13, G_23, TEMPERATURE

    const double E_1 = elastic_properties_.at(0);
    const double E_2 = elastic_properties_.at(1);
    const double nu_12 = elastic_properties_.at(2);
    const double G_12 = elastic_properties_.at(3);
    const double G_13 = elastic_properties_.at(4);
    const double G_23 = elastic_properties_.at(5);

    const double nu_21 = nu_12 * (E_2 / E_1);

    const double Q_11 = E_1 / (1.0 - nu_12 * nu_21);
    const double Q_12 = (nu_12 * E_2) / (1.0 - nu_12 * nu_21);
    const double Q_22 = E_2 / (1.0 - nu_12 * nu_21);
    const double Q_66 = G_12;
    const double Q_44 = G_23;
    const double Q_55 = G_13;

    std::vector<double> reduced_stiffness{Q_11, Q_12, Q_22, Q_66, Q_44, Q_55};

    return reduced_stiffness;
}