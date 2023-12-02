#ifndef FE_MODEL_H
#define FE_MODEL_H

#include <string>
#include <fstream>
#include <memory>
#include <map>

#include <Eigen/Sparse>

#include "Element.hpp"
#include "Material.hpp"
#include "Node.hpp"
#include "Section.hpp"
#include "Triplet.hpp"

/*
    TODO

    1) For better performance, use maps to store nodes and elements.
       - I'm not sure if this will indeed result in better performance.

    2) For improved performance, use an element's tag when storing element sets.
       - Solved (?): Using pointers to assign nodes to elements.

    3) Currently, the input file must follow a specific order of keywords, but
       the code should be able to read them in any order.
       - Is it a good idea to reopen the file multiple times to read a specific
         keyword?
*/

// Class that owns all model properties
class FEModel
{
public:
    FEModel(const std::string filename) : filename_(filename)
    {
        std::fstream input(filename_);
        std::string line;

        if (!input.is_open())
            throw std::runtime_error("\nINPUT FILE NOT FOUND!");
    };

    // The rule of five
    FEModel() = delete;
    FEModel(FEModel const &) = default;
    FEModel &operator=(FEModel const &) = default;
    FEModel(FEModel &&) = default;
    FEModel &operator=(FEModel &&) = default;

    //// Selectors
    inline std::size_t n_nodes() const { return node_map.size(); }
    inline std::size_t n_elements() const { return element_map.size(); }

    // TODO 
    // Considering 5 DOF per node here. If the same node is shared by
    // multiple types of elements, determine how to calculate the total number
    // of DOFs of the model.
    
    const std::size_t n_dof = 5;
    // Return the total number of DOF for the model.
    inline std::size_t total_dof() const { return n_dof * n_nodes(); }

    // Read mesh
    void read_input();

    void print_nodes();
    void print_elements();

private:
    // Smart pointers for defining containers of derived classes and shared
    // objects
    typedef std::shared_ptr<Boundary> BoundaryPtr;
    typedef std::shared_ptr<CLoad> CLoadPtr;
    typedef std::shared_ptr<Node> NodePtr;
    typedef std::shared_ptr<Material> MaterialPtr;
    typedef std::shared_ptr<Section> SectionPtr;
    typedef std::shared_ptr<DLoad> DLoadPtr;
    typedef std::shared_ptr<Element> ElementPtr;

public: // for debugging purposes
    const std::string filename_{};

    /*
        TO TEST AND STUDY

        1) How to store all mesh data? Just store the data itself and spread
           this data for each specific class after reading it or read it and
           defined all classes at once (currently)?

        2) Use vector, map or unordered map? Need some benchmarks, but unordered
           map should be the best option here.
    */

    // NODE TAG, BOUNDARY
    std::map<std::size_t, BoundaryPtr> boundary_map{};

    // NODE TAG, CLOAD
    std::map<std::size_t, CLoadPtr> cload_map{};

    // NODE TAG, NODE
    std::map<std::size_t, NodePtr> node_map{};

    // MATERIAL NAME, MATERIAL
    std::map<std::string, MaterialPtr> material_map{};

    // ELEMENT TAG, DLOAD
    std::map<std::size_t, DLoadPtr> dload_map{};

    // ELEMENT SET, SECTION
    std::map<std::string, SectionPtr> section_map{};

    // ELEMENT TAG, ELEMENT
    std::map<std::size_t, ElementPtr> element_map{};

    void read_boundary();
    void read_cload();
    void read_nodes();
    void read_materials();
    void read_sections();
    void read_dload();
    void read_elements();
    void read_temperature();

    // Check if the stiffness, or load, degree of freedom (dof) is in the
    // boundary condition dof.
    // Zero out all row and column entries that contain any boundary dof.
    // Set "X" to the diagonal entry.
    bool check_boundary_dof(const std::vector<std::size_t> dofs,
                            const std::size_t row) const;

    bool check_boundary_dof(const std::vector<std::size_t> dofs,
                            const std::size_t row,
                            const std::size_t col) const;

    //// Assembling

    // Stiffness
    std::vector<Triplet> 
    linear_stiffness_matrix() const;

    std::vector<Triplet> 
    nonlinear_stiffness_matrix(const std::vector<double> &solution) const;

    std::vector<Triplet>
    tangent_stiffness_matrix(const std::vector<double> &solution) const;

    // Force

    // Element distributed loading vector
    std::vector<Triplet> element_pressure_load_vector() const;

    // Element thermal loading vector
    std::vector<Triplet> element_thermal_load_vector() const;

    // Nodal force vector
    std::vector<Triplet> nodal_load_vector() const;

    // Total force vector
    std::vector<Triplet> force_vector() const;

    // Internal force vector
    // F_int = K(U)U
    std::vector<Triplet> 
    internal_force_vector(const std::vector<double> &solution) const;

    // Residual vector
    // R = K(U)U - F = F_int - F_ext
    std::vector<Triplet> 
    residual_vector(const std::vector<double> &solution) const;

    //// Solver
    std::vector<double> linear_solver() const;
    std::vector<double> nonlinear_solver() const;

    //// Output
    void write_result(const std::vector<double> &solution,
                      std::string filename) const;
};

#endif // FE_MODEL_H
