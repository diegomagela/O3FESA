#include "FEModel.hpp"
#include "FEModelIO.hpp"
#include "Elements.hpp"
#include "Materials.hpp"
#include "IO.hpp"
#include "Section.hpp"

#include <algorithm>
#include <fstream>

// TODO
//
// 1) Ignore comment lines in input file
// 2) Overall improvement: Be sure that keywords match exactly

FEModel::FEModel(const std::string filename)
{
    // Input stream
    std::fstream input(filename);

    // String to receive each files's line when reading it
    std::string line;

    if (input.is_open())
    {
        std::cout << "Reading the input file..." << '\n';

        while (true)
        {
            // Read the file's first line when opening the file
            // Read the next line otherwise
            std::getline(input, line);

            // Read nodes
            // Search for keyword "*NODE"
            // find() returns a size_t type, so a static_cast is used to compare
            if (line.find("*NODE") == static_cast<std::string::size_type>(0))
            {
                // Read line until reach the next section starting with '*'
                // or reach the end of the file
                //
                // peek() function "previews" the next character in the input
                while (input.peek() != '*' and !input.eof())
                {
                    std::getline(input, line);
                    add_node(line);
                }
            }

            // Read elements
            // Search for keyword "*ELEMENT"

            if (line.find("*ELEMENT") == static_cast<std::string::size_type>(0))
            {
                std::string element_set = get_element_set(line);
                std::string element_type = get_element_type(line);

                std::vector<std::size_t> element_set_vec{};

                while (input.peek() != '*' and !input.eof())
                {
                    std::getline(input, line);
                    add_element(line, element_type, element_set_vec);
                }

                element_sets.emplace(element_set, element_set_vec);
            }

            // Read boundary
            // Search for keyword "*BOUNDARY"

            if (line.find("*BOUNDARY") == static_cast<std::string::size_type>(0))
            {
                while (input.peek() != '*' and !input.eof())
                {
                    std::getline(input, line);
                    add_boundary(line);
                }
            }

            // Read cload (nodal loading)
            // Search for keyword "*CLOAD"

            if (line.find("*CLOAD") == static_cast<std::string::size_type>(0))
            {
                while (input.peek() != '*' and !input.eof())
                {
                    std::getline(input, line);
                    add_cload(line);
                }
            }

            // Read dload (distributed loading)
            // Search for keyword "*DLOAD"

            if (line.find("*DLOAD") == static_cast<std::string::size_type>(0))
            {
                while (input.peek() != '*' and !input.eof())
                {
                    std::getline(input, line);
                    add_dload(line);
                }
            }

            // Read material
            // Search for keyword "*MATERIAL"
            // TODO
            // 1) Check if the next line keyword is *ELASTIC
            // 2) Them check if *DENSITY is defined (although is not mandatory)

            if (line.find("*MATERIAL") == static_cast<std::string::size_type>(0))
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
                std::vector<double> elastic_properties_vec;

                for (const auto &properties : elastic_properties_str)
                    elastic_properties_vec.push_back(std::stod(properties));

                // Get density
                std::getline(input, line);
                std::getline(input, line);
                double density = std::stod(line);

                MaterialPtr material;

                if (material_type == "LAMINA")
                {
                    material = std::make_shared<Lamina>(material_name,
                                                        elastic_properties_vec,
                                                        density);
                }

                if (material_type == "ISOTROPIC")
                {
                    material = std::make_shared<Isotropic>(material_name,
                                                           elastic_properties_vec,
                                                           density);
                }

                material_map.emplace(material_name, material);
            }

            // Read shell section
            // Search for keyword "*SHELL SECTION"

            if (line.find("*SHELL SECTION") == static_cast<std::string::size_type>(0))
            {
                // By now, I will just consider COMPOSITE parameter
                std::string shell_section_set = get_shell_section_set(line);

                std::vector<double> thickness_vec{};
                std::vector<std::size_t> nip_vec{};
                std::vector<std::string> material_vec{};
                std::vector<int> orientation_vec{};

                double thickness{};
                std::size_t nip{};
                std::string material{};
                int orientation{};

                std::vector<std::string> shell_section_str{};

                while (input.peek() != '*' and !input.eof())
                {
                    std::getline(input, line);

                    shell_section_str = split_string(line, ',');

                    thickness = std::stod(shell_section_str.at(0));
                    nip = std::stoul(shell_section_str.at(1));
                    material = shell_section_str.at(2);
                    orientation = std::stoi(shell_section_str.at(3));

                    thickness_vec.push_back(thickness);
                    nip_vec.push_back(nip);
                    material_vec.push_back(material);
                    orientation_vec.push_back(orientation);
                }

                SectionPtr section = std::make_shared<Section>(thickness_vec,
                                                               nip_vec,
                                                               material_vec,
                                                               orientation_vec);

                section_map[shell_section_set] = section;

                set_section_to_element();
            }

            // Exit loop if reach end of file
            if (input.eof())
                break;
        }

        input.close();
    }

    else
        std::cout << "Unable to open the input file";
}

// Public members functions

void FEModel::print_nodes()
{
    std::cout << "*NODE" << std::endl;
    std::cout << "TAG \t X \t Y \t Z" << std::endl;


    for (const auto &[tag, node] : node_ptr_map)
        node.get()->print();
}

void FEModel::print_elements()
{
    std::cout << "*ELEMENT" << '\n';

    for (const auto &[tag, element] : element_map)
        element.get()->print();
}

// void FEModel::print_element_sets()
// {
//     std::cout << "*ELSET" << '\n';

//     for (const auto &[tag, element] : element_sets)
//     {
//         std::cout << tag << '\n';

//         for (const auto &x : element)
//             std::cout << x << '\t';

//         std::cout << '\n';
//     }
// }

// void FEModel::print_boundary()
// {
//     std::cout << "*BOUNDARY" << '\n';

//     for (const auto &[node_tag, dof] : boundary_map)
//     {
//         std::cout << node_tag << '\t';

//         for (const auto &x : dof)
//             std::cout << x << '\t';

//         std::cout << '\n';
//     }
// }

// void FEModel::print_cload()
// {
//     std::cout << "*CLOAD" << '\n';

//     for (const auto &[node_tag, dof] : cload_map)
//     {
//         std::cout << node_tag << '\t';

//         for (const auto &x : dof)
//             std::cout << x << '\t';

//         std::cout << '\n';
//     }
// }

// void FEModel::print_dload()
// {
//     std::cout << "*DLOAD" << '\n';

//     for (const auto &[element_tag, value] : dload_map)
//     {
//         std::cout << element_tag << '\t' << value;
//         std::cout << '\n';
//     }
// }

// Private members functions

void FEModel::add_node(const std::string &input)
{
    std::vector<std::string> node_str = split_string(input, ',');

    size_t tag = std::stoul(node_str.at(0));
    double x = std::stod(node_str.at(1));
    double y = std::stod(node_str.at(2));
    double z = std::stod(node_str.at(3));

    NodePtr node = std::make_shared<Node>(tag, x, y, z);
    node_ptr_map.emplace(tag, node);
}

void FEModel::add_element(const std::string &input,
                          const std::string &element_type,
                          std::vector<std::size_t> &element_set)
{
    std::vector<std::string> element_str = split_string(input, ',');

    std::size_t element_tag = std::stoul(element_str.at(0));

    element_set.push_back(element_tag);

    std::vector<NodePtr> element_nodes{};
    NodePtr node_ptr{};

    // Looping starting from the second element, as the first one is the
    // element tag
    for (std::size_t i = 1; i < element_str.size(); i++)
    {
        std::string node_tag_str = element_str.at(i);
        std::size_t node_tag = std::stoul(node_tag_str);

        node_ptr = node_ptr_map[node_tag];

        element_nodes.push_back(node_ptr);
    }

    ElementPtr element{};

    if (element_type == "S9")
    {
        element = std::make_unique<S9>(element_tag, element_nodes);
        element_map.emplace(element_tag, std::move(element));
    }

    if (element_type == "S4")
    {
        element = std::make_unique<S4>(element_tag, element_nodes);
        element_map.emplace(element_tag, std::move(element));
    }
}

void FEModel::add_boundary(const std::string &input)
{
    std::vector<bool> imposed_dofs(6);
    std::vector<double> imposed_values(6);

    std::vector<std::string> boundary_str = split_string(input, ',');
    std::size_t node_tag = std::stoul(boundary_str.at(0));
    std::string named_constraint = boundary_str.at(1);

    Boundary boundary{}; // empty declaration to fill map outside conditional

    if (named_constraint.compare("ENCASTRE") == 0)
    {
        // For better readability, vector is filled with zeros.
        std::fill(imposed_values.begin(), imposed_values.end(), 0.0);
        std::fill(imposed_dofs.begin(), imposed_dofs.end(), true);

        Boundary _boundary(node_tag, imposed_dofs, imposed_values);
        boundary = _boundary;
    }

    /*
        It looks like that is not necessary to have a boundary map. Maybe just
        assign Boundary to its respective Node
    */

    node_ptr_map[node_tag].get()->set_boundary(boundary);
}

void FEModel::add_cload(const std::string &input)
{
    std::vector<bool> loading_dofs(6);
    std::vector<double> loading_values(6);

    std::vector<std::string> cload_str = split_string(input, ',');
    std::size_t node_tag = std::stoul(cload_str.at(0));
    std::size_t dof = std::stoul(cload_str.at(1));
    double value = std::stod(cload_str.at(2));

    loading_dofs.at(dof - 1) = true;
    loading_values.at(dof - 1) = value;

    CLoad cload(node_tag, loading_dofs, loading_values);

    node_ptr_map[node_tag].get()->set_cload(cload);
}

void FEModel::add_dload(const std::string &input)
{
    std::vector<std::string> cload_str = split_string(input, ',');
    std::size_t element_tag = std::stoul(cload_str.at(0));
    double value = std::stod(cload_str.at(1));

    DLoad dload(element_tag, true, value);

    element_map[element_tag].get()->set_dload(dload);
}

void FEModel::set_section_to_element()
{
    for (const auto &[elset, elements] : element_sets)
    {
        auto section = section_map[elset];

        for (const auto &element_tag : elements)
        {
            element_map[element_tag].get()->set_section(section);
        }
    }
}