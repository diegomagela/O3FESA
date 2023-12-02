#ifndef ISOTROPIC_HPP
#define ISOTROPIC_HPP

#include "Material.hpp"

// Isotropic material class
class Isotropic : public Material
{
public:
    Isotropic(const std::string name,
              const std::vector<double> elastic_properties,
              const double density,
              const std::vector<double> thermal_expansion) : Material(name,
                                                                      type_,
                                                                      elastic_properties,
                                                                      density,
                                                                      thermal_expansion),
                                                             elastic_properties_(elastic_properties),
                                                             thermal_expansion_(thermal_expansion){};
    ~Isotropic(){};

    // Selectors
    std::vector<double> elastic_coefficients() const override;

private:
    static constexpr inline MaterialType type_ = MaterialType::Isotropic;
    std::vector<double> elastic_properties_{};
    std::vector<double> thermal_expansion_{};
};

#endif // ISOTROPIC_HPP
