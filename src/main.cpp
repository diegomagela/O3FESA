#include <iostream>
#include <string>

#include "FEModel.hpp"

int main()
{
    FEModel model("mesh.inp");
    model.read_input();
    model.print_nodes();
    model.print_elements();

    return EXIT_SUCCESS;
}