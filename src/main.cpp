#include <iostream>
#include <string>

#include "FEModel.hpp"

int main()
{
    FEModel model("mesh.inp");

    model.print_nodes();
    // model.print_elements();
    // model.print_element_sets();
    // model.print_boundary();
    // model.print_cload();
    // model.print_dload();

    return EXIT_SUCCESS;
}