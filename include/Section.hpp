#ifndef SECTION
#define SECTION

#include "Material.hpp"

#include <memory>
#include <vector>
#include <string>

class Section
{
private:
    typedef std::shared_ptr<Material> MaterialPtr;

public:
    Section(const std::vector<double> thickness,
            const std::vector<std::size_t> nip,
            const std::vector<MaterialPtr> material,
            const std::vector<int> orientation) : thickness_(thickness),
                                                  nip_(nip),
                                                  material_(material),
                                                  orientation_(orientation){};

    // The rule of five
    Section() = default;
    Section(Section const &) = default;
    Section &operator=(Section const &) = default;
    Section(Section &&) = default;
    Section &operator=(Section &&) = default;

    // Modifiers
    inline void set_temperature(const double temperature) { temperature_ = temperature; }

    // Getter
    inline std::vector<MaterialPtr> get_materials() const { return material_; }

    // Selectors
    inline double total_thickness() const;
    inline std::size_t number_of_layers() const { return thickness_.size(); };
    std::vector<double> layer_position() const;
    std::vector<double> extensional_stiffness() const;
    std::vector<double> bending_stiffness() const;
    std::vector<double> bending_extensional_stiffness() const;
    std::vector<double> extensional_shear_stiffness() const;
    std::vector<double> thermal_force_resultants() const;
    std::vector<double> thermal_moment_resultants() const;

    // Return if the temperature is greater than 0
    bool has_temperature() const;

private:
    std::vector<double> transformed_reduced_stiffness(double orientation,
                                                      std::vector<double> elastic_coefficients) const;

    std::vector<double> transformed_thermal_coefficients(double orientation,
                                                         std::vector<double> thermal_coefficients) const;

private:
    std::vector<double> thickness_{};     // Thickness of each layer
    std::vector<std::size_t> nip_{};      // Number of integration points
    std::vector<MaterialPtr> material_{}; // List of materials
    std::vector<int> orientation_{};      // List of angle of each layer
    double temperature_{};                 // Temperature gradient
};

#endif // SECTION
