#include "Isotropic.hpp"

std::vector<double> Isotropic::elastic_coefficients() const
{
    // E, \nu, G
    const double young_modulus = elastic_properties_.at(0);
    const double poisson = elastic_properties_.at(1);
    const double shear_modulus = (1.0 / 2.0) * young_modulus / (1 + poisson);

    const double delta = young_modulus / (1.0 - poisson * poisson);

    // Elastic coefficients
    const double Q_11 = delta;
    const double Q_12 = delta * poisson;
    const double Q_22 = Q_11;
    const double Q_66 = delta * shear_modulus;
    const double Q_44 = Q_66;
    const double Q_55 = Q_66;

    std::vector<double> reduced_stiffness{Q_11, Q_12, Q_22, Q_66, Q_44, Q_55};

    return reduced_stiffness;
}