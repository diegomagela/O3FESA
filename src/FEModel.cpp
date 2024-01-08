#include "FEModel.hpp"
#include "FEModelIO.hpp"
#include "Elements.hpp"
#include "Materials.hpp"
#include "IO.hpp"
#include "VectorUtil.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>

void FEModel::read_input()
{
    const auto start{std::chrono::high_resolution_clock::now()};
    std::cout << "[Reading input file";

    // read_step();

    read_boundary();
    read_cload();
    read_nodes();
    read_materials();
    read_sections();
    read_temperature();
    read_dload();
    read_elements();

    const auto end{std::chrono::high_resolution_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};

    std::cout << " - " << elapsed_seconds << "]" << std::endl;
}

void FEModel::read_boundary()
{
    // Input stream
    std::fstream input(filename_);

    // String to receive each files's line when reading it
    std::string line{};

    bool find{false};

    while (!find)
    {
        std::getline(input, line);

        if (find_keyword(line, "*BOUNDARY"))
        {
            // Read line until reach the next section starting with '*'
            // or reach the end of the file
            // peek() function "previews" the next character in the input

            while (input.peek() != '*' and !input.eof())
            {
                std::getline(input, line);

                std::vector<bool> imposed_dofs(6);
                std::vector<double> imposed_values(6);

                std::vector<std::string> boundary_str = split_string(line, ',');
                std::size_t node_tag = std::stoul(boundary_str.at(0));
                std::string named_constraint = boundary_str.at(1);

                if (!boundary_map.contains(node_tag))
                {
                    BoundaryPtr boundary{}; // empty declaration to fill map outside conditional

                    if (named_constraint.compare("ENCASTRE") == 0)
                    {
                        // For better readability, vector is filled with zeros.
                        std::fill(imposed_values.begin(), imposed_values.end(), 0.0);
                        std::fill(imposed_dofs.begin(), imposed_dofs.end(), true);

                        boundary = std::make_shared<Boundary>(node_tag, imposed_dofs, imposed_values);
                        boundary_map.emplace(node_tag, boundary);
                    }

                    else
                    {
                        std::size_t boundary_dof = std::stoul(boundary_str.at(1));
                        imposed_dofs.at(boundary_dof - 1) = true;

                        boundary = std::make_shared<Boundary>(node_tag, imposed_dofs, imposed_values);
                        boundary_map.emplace(node_tag, boundary);
                    }
                }

                else
                    boundary_map[node_tag].get()->set_imposed_dof(std::stoul(named_constraint) - 1);
            }
        }

        if (input.eof())
            find = true;
    }
}

void FEModel::read_cload()
{
    // Input stream
    std::fstream input(filename_);

    // String to receive each files's line when reading it
    std::string line{};

    bool find{false};

    while (!find)
    {
        std::getline(input, line);

        if (find_keyword(line, "*CLOAD"))
        {
            while (input.peek() != '*' and !input.eof())
            {
                std::getline(input, line);

                std::vector<bool> loading_dofs(6);
                std::vector<double> loading_values(6);

                std::vector<std::string> cload_str = split_string(line, ',');
                std::size_t node_tag = std::stoul(cload_str.at(0));
                std::size_t dof = std::stoul(cload_str.at(1));
                double value = std::stod(cload_str.at(2));

                loading_dofs.at(dof - 1) = true;
                loading_values.at(dof - 1) = value;

                CLoadPtr cload = std::make_shared<CLoad>(node_tag,
                                                         loading_dofs,
                                                         loading_values);

                cload_map.emplace(node_tag, cload);
            }
        }

        if (input.eof())
            find = true;
    }
}

void FEModel::read_nodes()
{
    // Input stream
    std::fstream input(filename_);

    // String to receive each files's line when reading it
    std::string line{};

    bool find{false};

    while (!find)
    {
        std::getline(input, line);

        if (find_keyword(line, "*NODE"))
        {
            // Read line until reach the next section starting with '*'
            // or reach the end of the file
            // peek() function "previews" the next character in the input

            while (input.peek() != '*' and !input.eof())
            {
                std::getline(input, line);

                std::vector<std::string> node_str = split_string(line, ',');

                size_t tag = std::stoul(node_str.at(0));
                double x = std::stod(node_str.at(1));
                double y = std::stod(node_str.at(2));
                double z = std::stod(node_str.at(3));

                NodePtr node = std::make_shared<Node>(tag, x, y, z);

                if (boundary_map.contains(tag))
                    node.get()->set_boundary(boundary_map[tag]);

                if (cload_map.contains(tag))
                    node.get()->set_cload(cload_map[tag]);

                node_map.emplace(tag, node);
            }
        }

        if (input.eof())
            find = true;
    }
}

void FEModel::read_materials()
{
    // Input stream
    std::fstream input(filename_);

    // String to receive each files's line when reading it
    std::string line{};

    bool find{false};

    while (!find)
    {
        std::getline(input, line);

        if (find_keyword(line, "*MATERIAL"))
        {
            // Get material name
            std::string material_name = get_material_name(line);

            // Get material type
            std::getline(input, line);
            std::string material_type = get_material_name(line);

            // Get elastic properties
            std::getline(input, line);

            std::vector<std::string> elastic_properties_str =
                split_string(line, ',');

            std::vector<double> elastic_properties_vec{};

            for (const auto &properties : elastic_properties_str)
                elastic_properties_vec.push_back(std::stod(properties));

            // Get density
            std::getline(input, line);
            std::getline(input, line);
            double density = std::stod(line);

            // Get thermal expansion coefficients
            std::vector<double> thermal_expansion_vec(2);
            std::getline(input, line);

            if (find_keyword(line, "*EXPANSION"))
            {
                std::string expansion_type = get_expansion_type(line);
                std::getline(input, line);

                if (expansion_type == "ORTHO")
                {
                    std::vector<std::string> thermal_expansion_str =
                        split_string(line, ',');

                    double alpha_1 = std::stod(thermal_expansion_str.at(0));
                    double alpha_2 = std::stod(thermal_expansion_str.at(1));

                    thermal_expansion_vec.at(0) = alpha_1;
                    thermal_expansion_vec.at(1) = alpha_2;
                }

                if (expansion_type == "ISO")
                {
                    std::vector<std::string> thermal_expansion_str =
                        split_string(line, ',');

                    if (thermal_expansion_str.size() > 1)
                    {
                        std::string const error =
                            "\nERROR: EXPECTED ONE PARAMETER FOR *EXPANSION, TYPE=ISO";

                        throw std::runtime_error(error);
                    }

                    double alpha = std::stod(thermal_expansion_str.at(0));

                    thermal_expansion_vec.at(0) = alpha;
                    thermal_expansion_vec.at(1) = alpha;
                }
            }

            MaterialPtr material;

            if (material_type == "LAMINA")
            {
                material = std::make_shared<Lamina>(material_name,
                                                    elastic_properties_vec,
                                                    density,
                                                    thermal_expansion_vec);
            }

            if (material_type == "ISOTROPIC")
            {
                material = std::make_shared<Isotropic>(material_name,
                                                       elastic_properties_vec,
                                                       density,
                                                       thermal_expansion_vec);
            }

            material_map.emplace(material_name, material);
        }

        if (input.eof())
            find = true;
    }
}

void FEModel::read_temperature()
{
    // Input stream
    std::fstream input(filename_);

    // String to receive each files's line when reading it
    std::string line{};

    bool find{false};

    while (!find)
    {
        std::getline(input, line);

        if (find_keyword(line, "*TEMPERATURE"))
        {
            while (input.peek() != '*' and !input.eof())
            {
                std::getline(input, line);

                std::vector<std::string> temperature_str = split_string(line, ',');
                std::string element_set = temperature_str.at(0);
                double temperature = std::stod(temperature_str.at(1));

                if (section_map.contains(element_set))
                    section_map[element_set].get()->set_temperature(temperature);

                else
                    throw std::runtime_error("\nERROR: ELEMENT SET DEFINED IN *TEMPERATURE NOT FOUND!");
            }
        }

        if (input.eof())
            find = true;
    }
}

void FEModel::read_sections()
{
    // Input stream
    std::fstream input(filename_);

    // String to receive each files's line when reading it
    std::string line{};

    bool find{false};

    while (!find)
    {
        std::getline(input, line);

        if (find_keyword(line, "*SHELL SECTION"))
        {
            // By now, I will just consider COMPOSITE parameter
            std::string shell_section_set = get_shell_section_set(line);

            std::vector<double> thickness_vec{};
            std::vector<std::size_t> nip_vec{};
            std::vector<MaterialPtr> material_vec{};
            std::vector<double> orientation_vec{};

            double thickness{};
            std::size_t nip{};
            std::string material_name{};
            int orientation{};

            std::vector<std::string> shell_section_str{};

            while (input.peek() != '*' and !input.eof())
            {
                std::getline(input, line);

                shell_section_str = split_string(line, ',');

                thickness = std::stod(shell_section_str.at(0));
                nip = std::stoul(shell_section_str.at(1));
                material_name = shell_section_str.at(2);
                orientation = std::stoi(shell_section_str.at(3));

                MaterialPtr material;

                if (material_map.contains(material_name))
                    material = material_map[material_name];

                else
                {
                    std::string error = "\n MATERIAL " +
                                        material_name +
                                        " DEFINED IN *SHELL SECTION NOT FOUND.\n";

                    throw std::runtime_error(error);
                }

                thickness_vec.push_back(thickness);
                nip_vec.push_back(nip);
                material_vec.push_back(material);
                orientation_vec.push_back(orientation);
            }

            SectionPtr section = std::make_shared<Section>(thickness_vec,
                                                           nip_vec,
                                                           material_vec,
                                                           orientation_vec);

            section_map.emplace(shell_section_set, section);
        }

        if (input.eof())
            find = true;
    }
}

void FEModel::read_dload()
{
    // Input stream
    std::fstream input(filename_);

    // String to receive each files's line when reading it
    std::string line{};

    bool find{false};

    while (!find)
    {
        std::getline(input, line);

        if (find_keyword(line, "*DLOAD"))
        {
            while (input.peek() != '*' and !input.eof())
            {
                std::getline(input, line);

                std::vector<std::string> cload_str = split_string(line, ',');
                std::size_t element_tag = std::stoul(cload_str.at(0));
                double value = std::stod(cload_str.at(1));

                DLoadPtr dload = std::make_shared<DLoad>(element_tag, value);

                dload_map.emplace(element_tag, dload);
            }
        }

        if (input.eof())
            find = true;
    }
}

void FEModel::read_elements()
{
    // Input stream
    std::fstream input(filename_);

    // String to receive each files's line when reading it
    std::string line{};

    bool find{false};

    while (!find)
    {
        std::getline(input, line);

        if (find_keyword(line, "*ELEMENT"))
        {
            std::string element_set = get_element_set(line);
            std::string element_type = get_element_type(line);

            while (input.peek() != '*' and !input.eof())
            {
                std::getline(input, line);

                std::vector<std::string> element_str = split_string(line, ',');

                std::size_t element_tag = std::stoul(element_str.at(0));

                std::vector<NodePtr> element_nodes{};

                // Looping starting from the second element, as the first one is the
                // element tag
                for (std::size_t i = 1; i < element_str.size(); ++i)
                {
                    std::string node_tag_str = element_str.at(i);
                    std::size_t node_tag = std::stoul(node_tag_str);

                    NodePtr node_ptr = node_map[node_tag];

                    element_nodes.push_back(node_ptr);
                }

                ElementPtr element{};

                if (element_type == "Q9")
                    element = std::make_shared<Q9>(element_tag);

                if (element_type == "Q4")
                    element = std::make_unique<Q4>(element_tag);

                // Set element's nodes
                element.get()->set_nodes(element_nodes);

                // Set element's dload (if it has it)
                if (dload_map.contains(element_tag))
                    element.get()->set_dload(dload_map[element_tag]);

                // Set element's section (mandatory)
                element.get()->set_section(section_map[element_set]);

                element_map.emplace(element_tag, element);
            }
        }

        if (input.eof())
            find = true;
    }
}

void FEModel::read_step()
{
    // Input stream
    std::fstream input(filename_);

    // String to receive each files's line when reading it
    std::string line{};

    bool find{false};

    while (!find)
    {
        std::getline(input, line);

        if (find_keyword(line, "*STEP"))
        {
            std::string analysis_type = get_step_analysis_type(line);
            std::cout << analysis_type << std::endl;

            // Step step(analysis_type);
            // step_map.push_back(step);

            read_boundary();
        }

        if (find_keyword(line, "*END STEP"))
        {
            std::cout << "END STEP" << std::endl;
            find = true;
        }
    }
}

std::vector<Triplet> FEModel::nodal_load_vector() const
{
    // Number of nodes containing loading
    std::size_t n_nodes = cload_map.size();

    // Maximum number of DOF
    std::size_t n_dof = 5;

    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.5 * n_nodes * n_dof;

    std::vector<Triplet> nodal_load;
    nodal_load.reserve(estimation_nnz);

    for (const auto &[tag, cload] : cload_map)
    {
        std::vector<std::size_t> global_dofs = cload.get()->global_dofs();
        std::vector<double> load_values = cload.get()->load_values();

        for (std::size_t i = 0; i < global_dofs.size(); ++i)
        {
            std::size_t dof = global_dofs[i];
            double value = load_values[i];
            Triplet t(dof, 0, value);
            nodal_load.push_back(t);
        }
    }

    return nodal_load;
}

std::vector<Triplet> FEModel::external_load_vector() const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> external_load{};
    external_load.reserve(estimation_nnz);

    for (const auto &[tag, element] : element_map)
    {
        std::vector<Triplet> f_ext = element.get()->external_load_triplet();

        external_load.insert(external_load.end(),
                             f_ext.begin(),
                             f_ext.end());
    }

    std::vector<Triplet> f_nodal = nodal_load_vector();

    external_load.insert(external_load.end(),
                         f_nodal.begin(),
                         f_nodal.end());

    return external_load;
}

std::vector<Triplet> FEModel::linear_stiffness_matrix() const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> stiffness_matrix{};
    stiffness_matrix.reserve(estimation_nnz);

    for (const auto &[tag, element] : element_map)
    {
        std::vector<Triplet> k_linear =
            element.get()->linear_stiffness_triplet();

        stiffness_matrix.insert(stiffness_matrix.end(),
                                k_linear.begin(),
                                k_linear.end());
    }

    return stiffness_matrix;
}

std::vector<Triplet> FEModel::nonlinear_stiffness_matrix() const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> stiffness_matrix{};
    stiffness_matrix.reserve(estimation_nnz);

    for (const auto &[tag, element] : element_map)
    {
        std::vector<Triplet> k_nonlinear =
            element.get()->nonlinear_stiffness_triplet();

        stiffness_matrix.insert(stiffness_matrix.end(),
                                k_nonlinear.begin(),
                                k_nonlinear.end());
    }

    return stiffness_matrix;
}

std::vector<Triplet> FEModel::internal_force_vector() const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> internal_force;
    internal_force.reserve(estimation_nnz);

    for (const auto &[tag, element] : element_map)
    {
        std::vector<Triplet> f_int = element.get()->internal_load_triplet();

        internal_force.insert(internal_force.end(),
                              f_int.begin(),
                              f_int.end());
    }

    return internal_force;
}

std::vector<Triplet> FEModel::tangent_stiffness_matrix() const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> stiffness_matrix{};
    stiffness_matrix.reserve(estimation_nnz);

    for (const auto &[tag, element] : element_map)
    {
        std::vector<Triplet> k_tangent =
            element.get()->tangent_stiffness_triplet();

        stiffness_matrix.insert(stiffness_matrix.end(),
                                k_tangent.begin(),
                                k_tangent.end());
    }

    return stiffness_matrix;
}

std::vector<double> FEModel::u_displacements() const
{
    std::vector<double> u_displacements;
    u_displacements.reserve(n_nodes());

    for (auto const &[tag, node] : node_map)
        u_displacements.push_back(node.get()->get_u_displacements());

    return u_displacements;
}

std::vector<double> FEModel::v_displacements() const
{
    std::vector<double> v_displacements;
    v_displacements.reserve(n_nodes());

    for (auto const &[tag, node] : node_map)
        v_displacements.push_back(node.get()->get_v_displacements());

    return v_displacements;
}

std::vector<double> FEModel::w_displacements() const
{
    std::vector<double> w_displacements;
    w_displacements.reserve(n_nodes());

    for (auto const &[tag, node] : node_map)
        w_displacements.push_back(node.get()->get_w_displacements());

    return w_displacements;
}

std::vector<double> FEModel::phi_x_displacements() const
{
    std::vector<double> phi_x_displacements;
    phi_x_displacements.reserve(n_nodes());

    for (auto const &[tag, node] : node_map)
        phi_x_displacements.push_back(node.get()->get_phi_x_displacements());

    return phi_x_displacements;
}

std::vector<double> FEModel::phi_y_displacements() const
{
    std::vector<double> phi_y_displacements;
    phi_y_displacements.reserve(n_nodes());

    for (auto const &[tag, node] : node_map)
        phi_y_displacements.push_back(node.get()->get_phi_y_displacements());

    return phi_y_displacements;
}

std::vector<double> FEModel::total_displacements() const
{
    std::vector<double> displacements{};
    displacements.reserve(total_dof());

    for (auto const &[tag, node] : node_map)
    {
        std::vector<double> node_displacements =
            node.get()->get_displacements();

        displacements.insert(displacements.end(),
                             node_displacements.begin(),
                             node_displacements.end());
    }

    return displacements;
}

void FEModel::update_displacement(const std::vector<double> &displacement_solution) const
{
    // Transfer solution displacements to nodes
    for (auto const &[tag, node] : node_map)
    {
        const std::size_t n_dof = node.get()->n_dof();
        std::size_t first_index = (tag - 1) * n_dof;
        std::size_t last_index = first_index + (n_dof - 1);

        std::vector<double> node_displacements =
            vct::extract_elements(displacement_solution, first_index, last_index);

        Dofs node_dof(node_displacements);

        node.get()->set_dofs(node_dof);
    }
}

void FEModel::output(const std::string &filename) const
{
    std::ofstream file(filename);

    for (auto const &[tag, node] : node_map)
    {
        double x = node.get()->get_x();
        double y = node.get()->get_y();
        double z = node.get()->get_z();
        double u = node.get()->get_u_displacements();
        double v = node.get()->get_v_displacements();
        double w = node.get()->get_w_displacements();
        double phi_x = node.get()->get_phi_x_displacements();
        double phi_y = node.get()->get_phi_y_displacements();

        file << x << ','
             << y << ','
             << z << ','
             << u << ','
             << v << ','
             << w << ','
             << phi_x << ','
             << phi_y << '\n';
    }

    file.close();
}

Eigen::SparseMatrix<double>
FEModel::triplet_to_sparse_vector(const std::vector<Triplet> &triplet) const
{
    const std::size_t size = total_dof();
    Eigen::SparseMatrix<double> sparse_vector(size, 1);
    sparse_vector.setFromTriplets(triplet.begin(), triplet.end());

    return sparse_vector;
}

Eigen::SparseMatrix<double>
FEModel::triplet_to_sparse_matrix(const std::vector<Triplet> &triplet) const
{
    const std::size_t size = total_dof();
    Eigen::SparseMatrix<double> sparse_matrix(size, size);
    sparse_matrix.setFromTriplets(triplet.begin(), triplet.end());

    return sparse_matrix;
}

std::vector<double>
FEModel::solver_cg(const Eigen::SparseMatrix<double> &A,
                   const Eigen::SparseMatrix<double> &b) const
{
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> solver;

    Eigen::VectorXd X = solver.compute(A).solve(b);

    return vct::eigen_vector_to_std(X);
}

std::vector<double>
FEModel::solver_lu(const Eigen::SparseMatrix<double> &A,
                   const Eigen::SparseMatrix<double> &b) const
{
    Eigen::SparseLU<Eigen::SparseMatrix<double>,
                    Eigen::COLAMDOrdering<int>>
        solver;

    solver.analyzePattern(A);
    solver.factorize(A);

    Eigen::VectorXd X = solver.solve(b);

    return vct::eigen_vector_to_std(X);
}

std::vector<double> FEModel::solver_qr(const Eigen::SparseMatrix<double> &A,
                                       const Eigen::SparseMatrix<double> &b) const
{
    Eigen::SparseQR<Eigen::SparseMatrix<double>,
                    Eigen::COLAMDOrdering<int>>
        solver;

    solver.analyzePattern(A);
    solver.factorize(A);

    Eigen::VectorXd X = solver.solve(b);

    return vct::eigen_vector_to_std(X);
}

std::vector<double>
FEModel::solver_bicgstab(const Eigen::SparseMatrix<double> &A,
                         const Eigen::SparseMatrix<double> &b) const
{
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

    solver.compute(A);
    Eigen::VectorXd X = solver.solve(b);

    std::cout << "Iterations:      " << solver.iterations() << std::endl;
    std::cout << "Estimated error: " << solver.error() << std::endl;

    X = solver.solve(b);

    return vct::eigen_vector_to_std(X);
}

void FEModel::linear_solver() const
{
    // Assemble stiffness matrix
    Eigen::SparseMatrix<double> K =
        triplet_to_sparse_matrix(linear_stiffness_matrix());

    // Assemble load vector
    Eigen::SparseMatrix<double> F =
        triplet_to_sparse_vector(external_load_vector());

    // >>>
    // Check matrix rank and if it is invertible

    auto stiffness = K.toDense();
    Eigen::FullPivLU<Eigen::MatrixXd> lu(stiffness);
    std::string is_invertible{"Not!"};
    std::string is_symmetric{"Not!"};
    std::string rank_deficient{"Not!"};

    if (lu.isInvertible())
        is_invertible = "Yes!";

    if (stiffness.isApprox(stiffness.transpose()))
        is_symmetric = "Yes!";

    std::cout << "["
              << "Symmetric: " << is_symmetric
              << " | "
              << "Invertible: " << is_invertible
              << " | "
              << "Dof: " << total_dof()
              << " | "
              << "Rank: " << lu.rank()
              << "]"
              << std::endl;

    // Solve linear system
    std::vector<double> x = solver_lu(K, F);
    update_displacement(x);

    // Output result linear system
    const std::string filename = "../python/solution/linear.txt";
    output(filename);
}

void FEModel::nonlinear_solver_newton_raphson() const
{
    std::ofstream load_displacement{"../python/solution/load_displacement.txt"};

    const std::size_t n_steps{10};
    const double load_increment{1.0 / n_steps};
    double lambda{load_increment};

    // Constant terms during steps

    // External load F_ext
    std::vector<Triplet> external_load = external_load_vector();
    Eigen::SparseVector<double> f_ext = triplet_to_sparse_vector(external_load);

    // Linear stiffness matrix
    std::vector<Triplet> linear_matrix = linear_stiffness_matrix();
    Eigen::SparseMatrix<double> k_linear = triplet_to_sparse_matrix(linear_matrix);

    for (std::size_t step = 0; step < n_steps; ++step)
    {
        std::cout << "["
                  << "Load factor: "
                  << lambda * 100 << "%"
                  << "]"
                  << std::endl;

        // Calculate residual
        std::vector<Triplet> internal_load = internal_force_vector();
        Eigen::SparseVector<double> f_int = triplet_to_sparse_vector(internal_load);
        Eigen::SparseVector<double> delta_f = lambda * f_ext;
        Eigen::SparseVector<double> residual = delta_f + f_int;
        double residual_norm = residual.norm();
        double residual_error = residual_norm /
                                std::max(delta_f.norm(), f_int.norm());

        // Calculate tangent matrix
        std::vector<Triplet> tangent_matrix = tangent_stiffness_matrix();
        Eigen::SparseMatrix<double> k_tangent = triplet_to_sparse_matrix(tangent_matrix);
        Eigen::SparseMatrix<double> k_total = k_tangent + k_linear;

        // Calculate displacement increment
        std::vector<double> dx = solver_lu(k_total, residual);

        // Update solution
        Eigen::VectorXd X0 = vct::std_to_eigen_vector(total_displacements());
        Eigen::VectorXd DX = vct::std_to_eigen_vector(dx);
        Eigen::VectorXd X1 = X0 + DX;
        update_displacement(vct::eigen_vector_to_std(X1));
        double displacement_error = (X1 - X0).norm() / X1.norm();

        std::size_t cont{1};
        std::size_t cont_tan{1};

        // Recalculate tangent matrix
        tangent_matrix = tangent_stiffness_matrix();
        k_tangent = triplet_to_sparse_matrix(tangent_matrix);
        k_total = k_linear + k_tangent;

        // >>>
        std::cout << "[Tangent norm: " << k_total.norm() << "]" << std::endl;

        // >>>
        // Check matrix rank and if it is invertible

        // auto stiffness = k_total.toDense();
        // Eigen::FullPivLU<Eigen::MatrixXd> lu(stiffness);
        // std::string is_invertible{"Not!"};
        // std::string is_symmetric{"Not!"};
        // std::string rank_deficient{"Not!"};

        // if (lu.isInvertible())
        //     is_invertible = "Yes!";

        // if (stiffness.isApprox(stiffness.transpose()))
        //     is_symmetric = "Yes!";

        // std::cout << "["
        //           << "Symmetric: " << is_symmetric
        //           << " | "
        //           << "Invertible: " << is_invertible
        //           << " | "
        //           << "Dof: " << total_dof()
        //           << " | "
        //           << "Rank: " << lu.rank()
        //           << "]"
        //           << std::endl;

        while (residual_error > 1E-3 or displacement_error > 1E-3)
        {
            // Recalculate tangent matrix
            // tangent_matrix = tangent_stiffness_matrix();
            // k_tangent = triplet_to_sparse_matrix(tangent_matrix);
            // k_total = k_tangent + k_linear;

            // Recalculate residual
            internal_load = internal_force_vector();
            f_int = triplet_to_sparse_vector(internal_load);

            std::cout << "[Internal force norm: " << f_int.norm() << "]" << std::endl;

            residual = delta_f + f_int;
            residual_norm = residual.norm();
            residual_error = residual_norm /
                             std::max(delta_f.norm(), f_int.norm());

            // Calculate increment
            dx = solver_lu(k_total, residual);

            // Update solution
            X0 = vct::std_to_eigen_vector(total_displacements());
            DX = vct::std_to_eigen_vector(dx);
            X1 = X0 + DX;
            update_displacement(vct::eigen_vector_to_std(X1));
            displacement_error = (X1 - X0).norm() / X1.norm();

            if (std::isnan(displacement_error) or std::isnan(residual_error))
            {
                std::cout << "[Solution diverges]" << std::endl;
                return;
            }

            std::cout << "["
                      << "Step: " << cont
                      << " | "
                      << "Displacement error: " << displacement_error
                      << " | "
                      << "Residual error: " << residual_error
                      << "]"
                      << std::endl;

            if (cont_tan > 3)
            {
                // Recalculate tangent matrix for better convergence
                tangent_matrix = tangent_stiffness_matrix();
                k_tangent = triplet_to_sparse_matrix(tangent_matrix);
                k_total = k_tangent + k_linear;

                std::cout << "["
                          << "Recalculating tangent matrix"
                          << "]"
                          << std::endl;

                cont_tan = 0;

                // >>>
                std::cout << "[Tangent norm: " << k_total.norm() << "]" << std::endl;
                std::cout << "[Internal force norm: " << f_int.norm() << "]" << std::endl;

                // >>>
                // Check matrix rank and if it is invertible

                // auto stiffness = k_total.toDense();
                // Eigen::FullPivLU<Eigen::MatrixXd> lu(stiffness);
                // std::string is_invertible{"Not!"};
                // std::string is_symmetric{"Not!"};
                // std::string rank_deficient{"Not!"};

                // if (lu.isInvertible())
                //     is_invertible = "Yes!";

                // if (stiffness.isApprox(stiffness.transpose()))
                //     is_symmetric = "Yes!";

                // std::cout << "["
                //           << "Symmetric: " << is_symmetric
                //           << " | "
                //           << "Invertible: " << is_invertible
                //           << " | "
                //           << "Dof: " << total_dof()
                //           << " | "
                //           << "Rank: " << lu.rank()
                //           << "]"
                //           << std::endl;
            }

            ++cont;
            ++cont_tan;
        }

        std::string filename = "../python/solution/nonlinear.txt";
        output(filename);

        std::cout << "[Solution converged]" << std::endl
                  << std::endl;

        auto w_total = w_displacements();
        auto W = vct::std_to_eigen_vector(w_total);

        load_displacement << lambda << ',' << W.cwiseAbs().maxCoeff() << '\n';

        lambda += load_increment;

        // wait_on_enter();
    }
}

void FEModel::nonlinear_solver_newton_raphson(const std::size_t n_iterations) const
{
    std::ofstream load_displacement{"../python/solution/load_displacement.txt"};

    const double load_increment{1.0 / n_iterations};
    double lambda{load_increment};

    // Constant terms during steps

    // External load F_ext
    std::vector<Triplet> external_load = external_load_vector();
    Eigen::SparseVector<double> f_ext = triplet_to_sparse_vector(external_load);

    // Linear stiffness matrix
    std::vector<Triplet> linear_matrix = linear_stiffness_matrix();
    Eigen::SparseMatrix<double> k_linear = triplet_to_sparse_matrix(linear_matrix);

    for (std::size_t step = 0; step < n_iterations; ++step)
    {
        std::cout << "["
                  << "Load factor: "
                  << lambda * 100 << "%"
                  << "]"
                  << std::endl;

        // Calculate residual
        std::vector<Triplet> internal_load = internal_force_vector();
        Eigen::SparseVector<double> f_int = triplet_to_sparse_vector(internal_load);
        Eigen::SparseVector<double> delta_f = lambda * f_ext;
        Eigen::SparseVector<double> residual = delta_f + f_int;
        double residual_norm = residual.norm();
        double residual_error = residual_norm /
                                std::max(delta_f.norm(), f_int.norm());

        // Calculate tangent matrix
        std::vector<Triplet> tangent_matrix = tangent_stiffness_matrix();
        Eigen::SparseMatrix<double> k_tangent = triplet_to_sparse_matrix(tangent_matrix);
        Eigen::SparseMatrix<double> k_total = k_tangent + k_linear;

        // Calculate displacement increment
        std::vector<double> dx = solver_lu(k_total, residual);

        // Update solution
        Eigen::VectorXd X0 = vct::std_to_eigen_vector(total_displacements());
        Eigen::VectorXd DX = vct::std_to_eigen_vector(dx);
        Eigen::VectorXd X1 = X0 + DX;
        update_displacement(vct::eigen_vector_to_std(X1));
        double displacement_error = (X1 - X0).norm() / X1.norm();

        std::size_t cont{1};
        std::size_t cont_tan{1};

        // Recalculate tangent matrix
        tangent_matrix = tangent_stiffness_matrix();
        k_tangent = triplet_to_sparse_matrix(tangent_matrix);
        k_total = k_linear + k_tangent;

        // >>>
        // Check matrix rank and if it is invertible

        // auto stiffness = k_total.toDense();
        // Eigen::FullPivLU<Eigen::MatrixXd> lu(stiffness);
        // std::string is_invertible{};

        // if (lu.isInvertible())
        //     is_invertible = "Yes!";
        // else
        //     is_invertible = "Not!";

        // std::cout << "["
        //           << "Is stiffness matrix invertible? " << is_invertible
        //           << " | "
        //           << "rank(K) = " << lu.rank()
        //           << "]"
        //           << std::endl;

        while (residual_error > 1e-3 or displacement_error > 1e-3)
        {
            // Recalculate tangent matrix
            // tangent_matrix = tangent_stiffness_matrix();
            // k_tangent = triplet_to_sparse_matrix(tangent_matrix);
            // k_total = k_tangent + k_linear;

            // Recalculate residual
            internal_load = internal_force_vector();
            f_int = triplet_to_sparse_vector(internal_load);
            residual = delta_f + f_int;
            residual_norm = residual.norm();
            residual_error = residual_norm /
                             std::max(delta_f.norm(), f_int.norm());

            // Calculate increment
            dx = solver_lu(k_total, residual);

            // Update solution
            X0 = vct::std_to_eigen_vector(total_displacements());
            DX = vct::std_to_eigen_vector(dx);
            X1 = X0 + DX;
            update_displacement(vct::eigen_vector_to_std(X1));
            displacement_error = (X1 - X0).norm() / X1.norm();

            std::string filename = "../python/solution/nonlinear_" +
                                   std::to_string(cont) + ".txt";
            output(filename);

            if (std::isnan(displacement_error) or std::isnan(residual_error))
            {
                std::cout << "[Solution diverges]" << std::endl;
                return;
            }

            std::cout << "["
                      << "Step: " << cont
                      << " | "
                      << "Displacement error: " << displacement_error
                      << " | "
                      << "Residual error: " << residual_error
                      << "]"
                      << std::endl;

            if (cont_tan > 10)
            {
                // Recalculate tangent matrix for better convergence
                tangent_matrix = tangent_stiffness_matrix();
                k_tangent = triplet_to_sparse_matrix(tangent_matrix);
                k_total = k_tangent + k_linear;

                std::cout << "["
                          << "Recalculating tangent matrix"
                          << "]"
                          << std::endl;

                cont_tan = 0;

                // >>>
                // Check matrix rank and if it is invertible

                // auto stiffness = k_total.toDense();
                // Eigen::FullPivLU<Eigen::MatrixXd> lu(stiffness);
                // std::string is_invertible{};

                // if (lu.isInvertible())
                //     is_invertible = "Yes!";
                // else
                //     is_invertible = "Not!";

                // std::cout << "["
                //           << "Is stiffness matrix invertible? " << is_invertible
                //           << " | "
                //           << "rank(K) = " << lu.rank()
                //           << "]"
                //           << std::endl;
            }

            ++cont;
            ++cont_tan;
        }

        wait_on_enter();

        std::string filename = "../python/solution/nonlinear_" +
                               std::to_string(cont + 1) + ".txt";
        output(filename);

        std::cout << "[Solution converged]" << std::endl
                  << std::endl;

        auto w_total = w_displacements();
        auto W = vct::std_to_eigen_vector(w_total);

        load_displacement << lambda << ',' << W.cwiseAbs().maxCoeff() << '\n';

        lambda += load_increment;

        // wait_on_enter();
    }
}

void FEModel::nonlinear_solver_arc_length() const
{
    std::ofstream load_displacement{"../python/solution/load_displacement.txt"};

    // Number of steps
    const std::size_t n_steps{30};

    // Initial load increment \delta \lambda
    double delta_lambda{1};

    // Initial displacement
    Eigen::VectorXd X0 = vct::std_to_eigen_vector(total_displacements());

    // External load
    std::vector<Triplet> external_load = external_load_vector();
    Eigen::SparseVector<double> f_ext = triplet_to_sparse_vector(external_load);

    // Linear stiffness matrix
    std::vector<Triplet> linear_matrix = linear_stiffness_matrix();
    Eigen::SparseMatrix<double> k_linear = triplet_to_sparse_matrix(linear_matrix);

    // Residual vector
    std::vector<Triplet> internal_load = internal_force_vector();
    Eigen::SparseVector<double> f_int = triplet_to_sparse_vector(internal_load);
    Eigen::SparseVector<double> residual = f_ext + f_int;

    // Tangent matrix
    std::vector<Triplet> tangent_matrix = tangent_stiffness_matrix();
    Eigen::SparseMatrix<double> k_tangent = triplet_to_sparse_matrix(tangent_matrix);
    Eigen::SparseMatrix<double> k_total = k_tangent + k_linear;

    // Solve for u_hat and u_bar
    Eigen::VectorXd delta_u_hat = vct::std_to_eigen_vector(solver_lu(k_total, f_ext));
    Eigen::VectorXd delta_u_bar = vct::std_to_eigen_vector(solver_lu(k_total, residual));

    // Compute solution increment and update the solution
    Eigen::VectorXd delta_u = delta_u_bar + delta_lambda * delta_u_hat;
    Eigen::VectorXd Delta_u = delta_u;
    Eigen::VectorXd X = delta_u + X0;
    update_displacement(vct::eigen_vector_to_std(X));

    // Update load increment
    double lambda = delta_lambda;

    // Compute the arc-length
    double Delta_s = delta_lambda * delta_u_hat.norm();

    // Iteration

    // Update external load
    f_ext *= lambda;

    // Solve for u_hat
    delta_u_hat = vct::std_to_eigen_vector(solver_lu(k_total, f_ext));

    // Compute incremental load parameter, \delta \lambda, from quadratic equation
    // a (dl)^2 + b (dl) + c = 0

    double a = delta_u_hat.dot(delta_u_hat);
    double b = 2.0 * (delta_u_bar + Delta_u).dot(delta_u_hat);
    double c = (delta_u_bar + Delta_u).dot(delta_u_bar + Delta_u) - Delta_s * Delta_s;

    double delta_lambda_1 = (-b + std::sqrt(b * b - 4 * a * c)) / (2 * a);
    double delta_lambda_2 = (-b - std::sqrt(b * b - 4 * a * c)) / (2 * a);
}