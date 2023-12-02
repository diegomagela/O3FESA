#ifndef LAMINA_HPP
#define LAMINA_HPP

#include "Material.hpp"

// Lamina (orthotropic) material class
class Lamina : public Material
{
public:
    Lamina(const std::string name,
           const std::vector<double> elastic_properties,
           const double density,
           const std::vector<double> thermal_expansion) : Material(name,
                                                                   type_,
                                                                   elastic_properties,
                                                                   density,
                                                                   thermal_expansion),
                                                          elastic_properties_(elastic_properties),
                                                          thermal_expansion_(thermal_expansion){};
    ~Lamina(){};

    // Selectors
    std::vector<double> elastic_coefficients() const override;

private:
    inline static constexpr MaterialType type_ = MaterialType::Lamina;
    std::vector<double> elastic_properties_{};
    std::vector<double> thermal_expansion_{};
};

#endif // LAMINA_HPP
