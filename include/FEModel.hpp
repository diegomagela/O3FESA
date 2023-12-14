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
#include "Time.hpp"
#include "Triplet.hpp"

/*
        TODO

        1) How to store all mesh data? Just store the data itself and spread
           this data for each specific class after reading it or read it and
           defined all classes at once (currently)?

        2) Use vector, map or unordered map? Need some benchmarks, but unordered
           map should be the best option here.

        3) Ignore comments in input file

        4) Be sure that keywords match exactly as they should

        Derived class must be a pointer:
        - Material
        - Section
        - Element

        Why consider other classes as pointers as well?
*/

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

    // Number of dofs per node
    // **It is being considered that all nodes have the same number of dofs when
    // defining the model's total number of dofs. If the model has nodes shared
    // by different kinds of elements, how to properly calculate the total
    // number of dofs?
    const std::size_t n_dof = 5;

    // Inline function

    // Returns the total number of nodes
    inline std::size_t n_nodes() const { return node_map.size(); }
    // Return the total number of elements
    inline std::size_t n_elements() const { return element_map.size(); }
    // Return the total number of model's DOF.
    inline std::size_t total_dof() const { return n_dof * n_nodes(); }

    // Read mesh input file
    void read_input();

public: // Public for debugging purposes, it must be private
    // Definitions

    typedef std::shared_ptr<Boundary> BoundaryPtr;
    typedef std::shared_ptr<CLoad> CLoadPtr;
    typedef std::shared_ptr<Node> NodePtr;
    typedef std::shared_ptr<Material> MaterialPtr;
    typedef std::shared_ptr<Section> SectionPtr;
    typedef std::shared_ptr<DLoad> DLoadPtr;
    typedef std::shared_ptr<Element> ElementPtr;

    // Maps for storing data

    // Node tag, Boundary object
    std::map<std::size_t, BoundaryPtr> boundary_map{};
    // Node tag, CLoad object
    std::map<std::size_t, CLoadPtr> cload_map{};
    // Node tag, Node object
    std::map<std::size_t, NodePtr> node_map{};
    // Material name, Material object
    std::map<std::string, MaterialPtr> material_map{};
    // Element tag, DLoad object
    std::map<std::size_t, DLoadPtr> dload_map{};
    // Element set tag, Section object
    std::map<std::string, SectionPtr> section_map{};
    // Element tag, Element object
    std::map<std::size_t, ElementPtr> element_map{};

    // Functions to read keywords from the input file

    // Read Boundary objects from input
    void read_boundary();
    // Read CLoad objects from input
    void read_cload();
    // Read Node objects from input
    void read_nodes();
    // Read Material objects from input
    void read_materials();
    // Read temperature from input
    void read_temperature();
    // Read Section objects from input
    void read_sections();
    // Read DLoad objects from input
    void read_dload();
    // Read Element objects from input
    void read_elements();

    // Selectors

    // Total nodal force vector
    std::vector<Triplet> nodal_load_vector() const;
    // Global external load vector
    std::vector<Triplet> external_load_vector() const;
    // Global linear stiffness matrix
    std::vector<Triplet> linear_stiffness_matrix() const;
    // Global nonlinear stiffness matrix
    std::vector<Triplet> nonlinear_stiffness_matrix() const;
    // Global internal load vector
    std::vector<Triplet> internal_force_vector() const;
    // Global tangent stiffness matrix
    std::vector<Triplet> tangent_stiffness_matrix() const;

    // Returns model's u displacement
    std::vector<double> u_displacements() const;
    // Returns model's v displacement
    std::vector<double> v_displacements() const;
    // Returns model's w displacement
    std::vector<double> w_displacements() const;
    // Returns model's phi_x displacement
    std::vector<double> phi_x_displacements() const;
    // Returns model's phi_y displacement
    std::vector<double> phi_y_displacements() const;


    // Returns the total displacements of all nodes and dofs
    std::vector<double> total_displacements() const;
    // Receive a solution displacement vector and distributed it to each node
    void update_displacement(const std::vector<double> &solution) const;
    // Write out the result to a file
    void output(const std::string &filename) const;
    // Convert std::vector<Triplet> to Eigen::SparseMatrix<double> vector
    Eigen::SparseMatrix<double>
    triplet_to_sparse_vector(const std::vector<Triplet> &triplet) const;
    // Convert std::vector<Triplet> to Eigen::SparseMatrix<double> matrix
    Eigen::SparseMatrix<double>
    triplet_to_sparse_matrix(const std::vector<Triplet> &triplet) const;
    // CG solver
    std::vector<double> solver_cg(const Eigen::SparseMatrix<double> &A,
                                  const Eigen::SparseMatrix<double> &b) const;
    // LU solver
    std::vector<double> solver_lu(const Eigen::SparseMatrix<double> &A,
                                  const Eigen::SparseMatrix<double> &b) const;

    // Solvers

    // Linear solver
    void linear_solver() const;
    // Nonlinear solvers
    // - Newton-Raphson
    void nonlinear_solver_newton_raphson() const;

    // >>>
    void wait_on_enter() const
    {
        std::string dummy;
        std::cout << "Enter to continue..." << std::endl;
        std::getline(std::cin, dummy);
    }

private:
    const std::string filename_{};
};

#endif // FE_MODEL_H
