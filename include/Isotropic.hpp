#ifndef ISOTROPIC_HPP
#define ISOTROPIC_HPP

#include "Material.hpp"

// Isotropic material class
class Isotropic : public Material
{
public:
    Isotropic(const std::string name,
              const std::vector<double> elastic_properties,
              const double density) : Material(name,
                                               type_,
                                               elastic_properties,
                                               density),
                                      elastic_properties_(elastic_properties){};
    ~Isotropic(){};

    // Selectors
    std::vector<double> elastic_coefficients() const override;

private:
    static constexpr inline MaterialType type_ = MaterialType::Isotropic;
    std::vector<double> elastic_properties_{};
};

#endif // ISOTROPIC_HPP
