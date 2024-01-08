#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <vector>
#include <string>

enum class MaterialType
{
    Isotropic,
    Lamina
};

// Base class for other materials
class Material
{
public:
    Material(const std::string name,
             const MaterialType type,
             const std::vector<double> elastic_properties,
             const double density,
             const std::vector<double> thermal_expansion_coefficients) : name_(name),
                                                                         type_(type),
                                                                         elastic_properties_(elastic_properties),
                                                                         density_(density),
                                                                         thermal_expansion_coefficients_(thermal_expansion_coefficients){};

    // The rule of five
    Material() = default;
    Material(Material const &) = default;
    Material &operator=(Material const &) = default;
    Material(Material &&) = default;
    Material &operator=(Material &&) = default;

    // Selectors
    virtual std::vector<double> elastic_coefficients() const = 0;
    inline double density() const { return density_; };
    inline std::vector<double> thermal_expansion_coefficients() const
    {
        return thermal_expansion_coefficients_;
    }

private:
    std::string name_{};
    MaterialType type_{};
    std::vector<double> elastic_properties_{};
    double density_{};
    std::vector<double> thermal_expansion_coefficients_{};
};

#endif // MATERIAL_HPP
