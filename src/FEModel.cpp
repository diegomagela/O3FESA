#include "FEModel.hpp"
#include "FEModelIO.hpp"
#include "Elements.hpp"
#include "Materials.hpp"
#include "IO.hpp"
#include "VectorUtil.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>

// TODO
//
// 1) Ignore comment lines in input file
// 2) Be sure that keywords match exactly

void FEModel::read_input()
{
    std::cout << "Reading input file...";

    read_boundary();
    read_cload();
    read_nodes();
    read_materials();
    read_sections();
    read_temperature();
    read_dload();
    read_elements();

    std::cout << " done!" << '\n';
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

                BoundaryPtr boundary; // empty declaration to fill map outside conditional

                if (named_constraint.compare("ENCASTRE") == 0)
                {
                    // For better readability, vector is filled with zeros.
                    std::fill(imposed_values.begin(), imposed_values.end(), 0.0);
                    std::fill(imposed_dofs.begin(), imposed_dofs.end(), true);

                    boundary = std::make_shared<Boundary>(node_tag, imposed_dofs, imposed_values);
                }

                boundary_map.emplace(node_tag, boundary);
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
            std::vector<int> orientation_vec{};

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
                NodePtr node_ptr{};

                // Looping starting from the second element, as the first one is the
                // element tag
                for (std::size_t i = 1; i < element_str.size(); ++i)
                {
                    std::string node_tag_str = element_str.at(i);
                    std::size_t node_tag = std::stoul(node_tag_str);

                    node_ptr = node_map[node_tag];

                    element_nodes.push_back(node_ptr);
                }

                ElementPtr element{};

                if (element_type == "Q9")
                    element = std::make_shared<Q9>(element_tag);

                // if (element_type == "Q4")
                //     element = std::make_unique<Q4>(element_tag);

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

bool FEModel::check_boundary_dof(const std::vector<std::size_t> &dofs,
                                 const std::size_t row) const
{
    if (vct::vector_has_element(dofs, row))
        return true;

    else
        return false;
}

bool FEModel::check_boundary_dof(const std::vector<std::size_t> &dofs,
                                 const std::size_t row,
                                 const std::size_t col) const
{
    if (vct::vector_has_element(dofs, row) or vct::vector_has_element(dofs, col))
        return true;

    else
        return false;
}

// =================================
// GLOBAL MODEL MATRICES AND VECTORS
// =================================

std::vector<Triplet> FEModel::linear_stiffness_matrix() const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> stiffness;
    stiffness.reserve(estimation_nnz);

    for (const auto &[tag, element] : element_map)
    {
        std::vector<std::size_t> element_global_dofs = element.get()->global_dofs();
        std::vector<std::size_t> element_local_dofs = element.get()->local_dofs();
        std::vector<std::size_t> boundary_dofs = element.get()->boundary_dofs();
        Eigen::MatrixXd element_stiffness = element.get()->linear_stiffness_matrix();

        for (std::size_t i = 0; i < element.get()->total_dof(); ++i)
        {
            for (std::size_t j = 0; j < element.get()->total_dof(); ++j)
            {
                // Local indexes
                std::size_t row_local = element_local_dofs.at(i);
                std::size_t col_local = element_local_dofs.at(j);

                // Global indexes
                std::size_t row_global = element_global_dofs.at(i);
                std::size_t col_global = element_global_dofs.at(j);

                double value = element_stiffness(row_local, col_local);

                if (check_boundary_dof(boundary_dofs, row_global, col_global))
                {
                    if (row_global == col_global)
                        value = element_stiffness.mean();

                    else
                        value = 0;
                }

                if (std::abs(value) > 0)
                {
                    Triplet t(row_global, col_global, value);
                    stiffness.push_back(t);
                }
            }
        }
    }

    return stiffness;
}

std::vector<Triplet> FEModel::nonlinear_stiffness_matrix(const std::vector<double> &solution) const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> stiffness;
    stiffness.reserve(estimation_nnz);

    for (const auto &[tag, element] : element_map)
    {
        std::vector<std::size_t> element_global_dofs = element.get()->global_dofs();
        std::vector<std::size_t> element_local_dofs = element.get()->local_dofs();
        std::vector<std::size_t> boundary_dofs = element.get()->boundary_dofs();
        std::vector<double> w_e = element.get()->get_w_displacements(solution);

        Eigen::MatrixXd element_stiffness =
            element.get()->nonlinear_stiffness_matrix(w_e) +
            element.get()->linear_stiffness_matrix();

        for (std::size_t i = 0; i < element.get()->total_dof(); ++i)
        {
            for (std::size_t j = 0; j < element.get()->total_dof(); ++j)
            {
                // Local indexes
                std::size_t row_local = element_local_dofs.at(i);
                std::size_t col_local = element_local_dofs.at(j);

                // Global indexes
                std::size_t row_global = element_global_dofs.at(i);
                std::size_t col_global = element_global_dofs.at(j);

                double value = element_stiffness(row_local, col_local);

                if (check_boundary_dof(boundary_dofs, row_global, col_global))
                {
                    if (row_global == col_global)
                        value = 1; // >>> CHANGE

                    else
                        value = 0;
                }

                if (std::abs(value) > 0)
                {
                    Triplet t(row_global, col_global, value);
                    stiffness.push_back(t);
                }
            }
        }
    }

    return stiffness;
}

std::vector<Triplet>
FEModel::tangent_stiffness_matrix(const std::vector<double> &solution) const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> stiffness;
    stiffness.reserve(estimation_nnz);

    for (const auto &[tag, element] : element_map)
    {
        std::vector<std::size_t> element_global_dofs = element.get()->global_dofs();
        std::vector<std::size_t> element_local_dofs = element.get()->local_dofs();
        std::vector<std::size_t> boundary_dofs = element.get()->boundary_dofs();
        std::vector<double> q_e = element.get()->get_displacements(solution);

        Eigen::MatrixXd
            element_stiffness = element.get()->tangent_stiffness_matrix(q_e);

        for (std::size_t i = 0; i < element.get()->total_dof(); ++i)
        {
            for (std::size_t j = 0; j < element.get()->total_dof(); ++j)
            {
                // Local indexes
                std::size_t row_local = element_local_dofs.at(i);
                std::size_t col_local = element_local_dofs.at(j);

                // Global indexes
                std::size_t row_global = element_global_dofs.at(i);
                std::size_t col_global = element_global_dofs.at(j);

                double value = element_stiffness(row_local, col_local);

                if (check_boundary_dof(boundary_dofs, row_global, col_global))
                {
                    if (row_global == col_global)
                        value = 1; // >>> CHANGE

                    else
                        value = 0;
                }

                if (std::abs(value) > 0)
                {
                    Triplet t(row_global, col_global, value);
                    stiffness.push_back(t);
                }
            }
        }
    }

    return stiffness;
}

std::vector<Triplet> FEModel::element_pressure_load_vector() const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> pressure_load;
    pressure_load.reserve(estimation_nnz);

    for (const auto &[tag, dload] : dload_map)
    {
        ElementPtr element = element_map.at(tag);

        std::vector<std::size_t> element_global_dofs = element.get()->global_dofs();
        std::vector<std::size_t> element_local_dofs = element.get()->local_dofs();
        std::vector<std::size_t> boundary_dofs = element.get()->boundary_dofs();
        Eigen::VectorXd load_vector = element.get()->pressure_load_vector();

        for (std::size_t i = 0; i < element.get()->total_dof(); ++i)
        {
            // Local index
            std::size_t row_local = element_local_dofs.at(i);

            // Global index
            std::size_t row_global = element_global_dofs.at(i);

            double value = load_vector(row_local);

            if (check_boundary_dof(boundary_dofs, row_global))
                value = 0;

            if (std::abs(value) > 0)
            {
                Triplet t(row_global, 0, value);
                pressure_load.push_back(t);
            }
        }
    }

    return pressure_load;
}

std::vector<Triplet> FEModel::element_thermal_load_vector() const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> thermal_load;
    thermal_load.reserve(estimation_nnz);

    for (const auto &[tag, element] : element_map)
    {
        if (element.get()->has_thermal_load())
        {

            std::vector<std::size_t> element_global_dofs = element.get()->global_dofs();
            std::vector<std::size_t> element_local_dofs = element.get()->local_dofs();
            std::vector<std::size_t> boundary_dofs = element.get()->boundary_dofs();
            Eigen::VectorXd load_vector = element.get()->thermal_load_vector();

            for (std::size_t i = 0; i < element.get()->total_dof(); ++i)
            {
                // Local index
                std::size_t row_local = element_local_dofs.at(i);

                // Global index
                std::size_t row_global = element_global_dofs.at(i);

                double value = load_vector(row_local);

                if (check_boundary_dof(boundary_dofs, row_global))
                    value = 0;

                if (std::abs(value) > 0)
                {
                    Triplet t(row_global, 0, value);
                    thermal_load.push_back(t);
                }
            }
        }
    }

    return thermal_load;
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

std::vector<Triplet> FEModel::force_vector() const
{
    std::vector<Triplet> element_pressure_load = element_pressure_load_vector();
    std::vector<Triplet> nodal_load = nodal_load_vector();
    std::vector<Triplet> element_thermal_load = element_thermal_load_vector();

    const std::size_t size = element_pressure_load.size() +
                             nodal_load.size() +
                             element_thermal_load.size();

    std::vector<Triplet> force_vector;
    force_vector.reserve(size);

    // Concatenate element and nodal loading vectors

    force_vector = nodal_load;

    force_vector.insert(force_vector.end(),
                        element_pressure_load.begin(),
                        element_pressure_load.end());

    force_vector.insert(force_vector.end(),
                        element_thermal_load.begin(),
                        element_thermal_load.end());

    return force_vector;
}

std::vector<Triplet> FEModel::internal_force_vector(const std::vector<double> &solution) const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.5 * total_dof();

    std::vector<Triplet> internal_force;
    internal_force.reserve(estimation_nnz);

    for (const auto &[tag, element] : element_map)
    {
        std::vector<std::size_t> element_global_dofs = element.get()->global_dofs();
        std::vector<std::size_t> element_local_dofs = element.get()->local_dofs();
        std::vector<std::size_t> boundary_dofs = element.get()->boundary_dofs();
        std::vector<double> d_e = element.get()->get_displacements(solution);

        Eigen::VectorXd element_internal_force = element.get()->internal_force(d_e);

        for (std::size_t i = 0; i < element.get()->total_dof(); ++i)
        {
            // Local index
            std::size_t row_local = element_local_dofs.at(i);

            // Global index
            std::size_t row_global = element_global_dofs.at(i);

            double value = element_internal_force(row_local);

            if (check_boundary_dof(boundary_dofs, row_global))
                value = 0;

            if (std::abs(value) > 0)
            {
                Triplet t(row_global, 0, value);
                internal_force.push_back(t);
            }
        }
    }

    return internal_force;
}

std::vector<Triplet>
FEModel::residual_vector(const std::vector<double> &solution) const
{
    std::vector<Triplet> F_ext = force_vector();
    std::vector<Triplet> F_int = internal_force_vector(solution);

    const std::size_t size = F_ext.size() + F_int.size();
    std::vector<Triplet> residual;
    residual.reserve(size);

    for (auto const &x : F_ext)
    {
        Triplet t(x.row(), x.col(), -1.0 * x.value());
        residual.push_back(t);
    }

    residual.insert(residual.end(), F_int.begin(), F_int.end());

    return residual;
}

std::vector<double> FEModel::linear_solver() const
{
    std::vector<Triplet> K = linear_stiffness_matrix();
    std::vector<Triplet> F = force_vector();

    Eigen::SparseMatrix<double> sm_K(total_dof(), total_dof());
    sm_K.setFromTriplets(K.begin(), K.end());

    Eigen::SparseMatrix<double> sm_F(total_dof(), 1);
    sm_F.setFromTriplets(F.begin(), F.end());

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper | Eigen::Lower> solver;
    Eigen::VectorXd X = solver.compute(sm_K).solve(sm_F);

    std::vector<double> solution;
    solution.resize(X.size());
    Eigen::VectorXd::Map(&solution[0], X.size()) = X;

    std::string filename = "../python/sol.txt";
    write_result(solution, filename);

    return solution;
}

std::vector<double> FEModel::nonlinear_solver() const
{
    // Initial solution from linear solution
    std::vector<double> x0 = linear_solver();

    // Fill T and F
    std::vector<Triplet> K = tangent_stiffness_matrix(x0);
    Eigen::SparseMatrix<double> sm_K(total_dof(), total_dof());
    sm_K.setFromTriplets(K.begin(), K.end());

    std::vector<Triplet> F = residual_vector(x0);
    Eigen::SparseMatrix<double> sm_F(total_dof(), 1);
    sm_F.setFromTriplets(F.begin(), F.end());

    // Solve
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(sm_K);
    Eigen::VectorXd X1 = solver.solve(sm_F);
    X1 = solver.solve(sm_F);

    Eigen::VectorXd X0 = Eigen::Map<Eigen::VectorXd>(x0.data(), x0.size());

    Eigen::VectorXd X = X0 + X1;

    // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    // solver.analyzePattern(sm_K);
    // solver.factorize(sm_K);
    // Eigen::VectorXd X = solver.solve(sm_F);

    // Output solution
    std::vector<double> solution;
    solution.resize(X.size());
    Eigen::VectorXd::Map(&solution[0], X.size()) = X;

    std::string filename = "/home/magela/Documents/Python/sol_nl_th.txt";
    write_result(solution, filename);

    return solution;
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

    for(auto const& [tag, node] : node_map)
    {
        double x = node.get()->get_x();
        double y = node.get()->get_y();
        double w = node.get()->get_w_displacements();

        file << x << ',' << y << ',' << w << '\n';
    }

    file.close();
}

void FEModel::write_result(const std::vector<double> &solution,
                           std::string filename) const
{
    const std::size_t size = n_nodes();

    std::vector<double> w_disp;
    w_disp.reserve(size);

    for (std::size_t i = 0; i < size; ++i)
        w_disp[i] = solution[5 * i + 2];

    std::ofstream file(filename);
    std::size_t cont = 0;
    for (auto const &[tag, node] : node_map)
    {
        double x = node.get()->get_x();
        double y = node.get()->get_y();
        double w = w_disp[cont];

        ++cont;

        file << x << ',' << y << ',' << w << '\n';
    }

    file.close();
}

// Public members functions

void FEModel::print_nodes()
{
    std::cout << "*NODE" << std::endl;

    for (const auto &[tag, node] : node_map)
        node.get()->print();
}

void FEModel::print_elements()
{
    std::cout << "*ELEMENT" << '\n';

    for (const auto &[tag, element] : element_map)
        element.get()->print();
}

//
//
//
//
//
//
//

std::vector<Triplet> FEModel::linear_stiffness_matrix_test() const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> stiffness;
    stiffness.reserve(estimation_nnz);

    for (const auto &[tag, element] : element_map)
    {
        vct::append_vector(stiffness,
                           element.get()->linear_stiffness_triplet());
    }

    return stiffness;
}

std::vector<Triplet> FEModel::element_load_vector_test() const
{
    // Estimated number of nonzero entries
    std::size_t estimation_nnz = 0.2 * total_dof();

    std::vector<Triplet> load;
    load.reserve(estimation_nnz);

    for (const auto &[tag, dload] : dload_map)
    {
        ElementPtr element = element_map.at(tag);
        vct::append_vector(load,
                           element.get()->element_load_triplet());
    }

    vct::append_vector(load, nodal_load_vector());

    return load;
}

std::vector<double> FEModel::linear_solver_test() const
{
    std::vector<Triplet> K = linear_stiffness_matrix_test();
    std::vector<Triplet> F = element_load_vector_test();

    Eigen::SparseMatrix<double> sm_K(total_dof(), total_dof());
    sm_K.setFromTriplets(K.begin(), K.end());

    Eigen::SparseMatrix<double> sm_F(total_dof(), 1);
    sm_F.setFromTriplets(F.begin(), F.end());

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper | Eigen::Lower> solver;
    Eigen::VectorXd X = solver.compute(sm_K).solve(sm_F);

    std::vector<double> solution;
    solution.resize(X.size());
    Eigen::VectorXd::Map(&solution[0], X.size()) = X;

    update_displacement(solution);
    output("test.txt");

    return solution;
}
