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
             double density) : name_(name),
                               type_(type),
                               elastic_properties_(elastic_properties),
                               density_(density){};

    // The rule of five
    Material() = default;
    Material(Material const &) = default;
    Material &operator=(Material const &) = default;
    Material(Material &&) = default;
    Material &operator=(Material &&) = default;

private:
    std::string name_{};
    MaterialType type_{};
    std::vector<double> elastic_properties_{};
    double density_{};
};

#endif // MATERIAL_HPP
