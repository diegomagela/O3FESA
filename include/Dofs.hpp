#ifndef DOFS_HPP
#define DOFS_HPP

#include <vector>

/*
        ==============
        DOF      INDEX
        ==============
        u        0
        v        1
        w        2
        phi_x    3
        phi_y    4
        phi_z    5
*/

// Class to handle nodes' displacements, velocities and accelerations
class Dofs
{
public:
    Dofs(std::vector<double> displacement) : displacements_(displacement) {}

    Dofs(std::vector<double> displacement,
         std::vector<double> velocities,
         std::vector<double> accelerations) : displacements_(displacement),
                                              velocities_(velocities),
                                              accelerations_(accelerations) {}

    // The rule of five

    Dofs() = default;
    Dofs(Dofs const &) = default;
    Dofs &operator=(Dofs const &) = default;
    Dofs(Dofs &&) = default;
    Dofs &operator=(Dofs &&) = default;

    // Selectors

    // Return dofs displacements
    inline std::vector<double>
    get_displacements() const { return displacements_; }
    // Return dofs velocities
    inline std::vector<double>
    get_velocities() const { return velocities_; }
    // Return dofs accelerations
    inline std::vector<double>
    get_accelerations() const { return accelerations_; }
    // Return u displacement
    inline double get_u_displacement() const { return displacements_.at(0); }
    // Return v displacement
    inline double get_v_displacement() const { return displacements_.at(1); }
    // Return w displacement
    inline double get_w_displacement() const { return displacements_.at(2); }
    // Return v displacement
    inline double get_phi_x_displacement() const { return displacements_.at(3); }
    // Return w displacement
    inline double get_phi_y_displacement() const { return displacements_.at(4); }

private:
    std::vector<double> displacements_{0, 0, 0, 0, 0};
    std::vector<double> velocities_{0, 0, 0, 0, 0};
    std::vector<double> accelerations_{0, 0, 0, 0, 0};
};

#endif // DOFS_HPP