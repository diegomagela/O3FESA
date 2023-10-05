#ifndef LAMINA_HPP
#define LAMINA_HPP

#include "Material.hpp"

// Lamina (orthotropic) material class
class Lamina : public Material
{
public:
    Lamina(const std::string name,
           const std::vector<double> elastic_properties,
           const double density) : Material(name,
                                            type_,
                                            elastic_properties,
                                            density){};
    ~Lamina(){};

    // Selectors
    std::vector<double> elastic_coefficients() const override;

private:
    inline static constexpr MaterialType type_ = MaterialType::Lamina;
    std::vector<double> elastic_properties_{};
};

#endif // LAMINA_HPP
