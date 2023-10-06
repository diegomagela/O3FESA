#include <iostream>
#include <string>

#include "FEModel.hpp"

int main()
{
    FEModel model("mesh.inp");
    model.read_input();
    model.print_nodes();
    model.print_elements();

    model.section_map["E1"].get()->extensional_stiffness();

    return EXIT_SUCCESS;
}