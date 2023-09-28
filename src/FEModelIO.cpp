#include "FEModelIO.hpp"
#include "IO.hpp"

std::string get_element_set(const std::string &input)
{
    std::vector<std::string> splitted_string = split_string(input, ',');

    /*
        TO-DO

     1) Implement error checker;

     2) Not necessarily the element type is at position 2, implement a search to
        find the correct position in the vector
    */

    // Check if an ELSET is defined, else define set as EALL
    std::string set{};

    if (splitted_string.size() < 3)
        set = "EALL";

    else
    {
        std::string element_set = splitted_string.at(2);
        set = split_string(element_set, '=').at(1);
    }

    return set;
}

std::string get_element_type(const std::string &input)
{
    std::vector<std::string> splitted_string = split_string(input, ',');
    std::string element_type = splitted_string.at(1);
    std::string type = split_string(element_type, '=').at(1);

    return type;
}

std::string get_material_name(const std::string &input)
{
    std::vector<std::string> splitted_string = split_string(input, ',');
    std::string material_name = splitted_string.at(1);
    std::string name = split_string(material_name, '=').at(1);

    return name;
}

std::string get_material_type(const std::string &input)
{
    std::vector<std::string> splitted_string = split_string(input, ',');
    std::string material_type = splitted_string.at(1);
    std::string type = split_string(material_type, '=').at(1);

    return type;
}

std::string get_shell_section_set(const std::string &input)
{
    std::vector<std::string> splitted_string = split_string(input, ',');
    std::string shell_section_set = splitted_string.at(1);
    std::string set = split_string(shell_section_set, '=').at(1);

    return set;
}