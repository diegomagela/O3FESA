#include "FEModel.hpp"
#include "FEModelIO.hpp"
#include "Elements.hpp"
#include "Materials.hpp"
#include "IO.hpp"

#include <algorithm>
#include <fstream>

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
                    std::cerr << '\n'
                              << "Material "
                              << material_name
                              << " defined in *SHELL SECTION not found!"
                              << '\n';
                              
                    throw std::exception();
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
                for (std::size_t i = 1; i < element_str.size(); i++)
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

// Public members functions

void FEModel::print_nodes()
{
    std::cout << "NODES" << std::endl;

    for (const auto &[tag, node] : node_map)
        node.get()->print();
}

void FEModel::print_elements()
{
    std::cout << "*ELEMENT" << '\n';

    for (const auto &[tag, element] : element_map)
        element.get()->print();
}