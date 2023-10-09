#include <iostream>
#include <string>

#include "FEModel.hpp"

template <class T>
void print_vector(std::vector<T> vector)
{
    for (auto const &x : vector)
        std::cout << x << '\t';

    std::cout << std::endl;
}

int main()
{
    FEModel model("mesh.inp");
    model.read_input();

    auto A = model.section_map["EALL"].get()->extensional_stiffness();
    auto B = model.section_map["EALL"].get()->bending_extensional_stiffness();
    auto D = model.section_map["EALL"].get()->bending_stiffness();
    auto C = model.section_map["EALL"].get()->extensional_shear_stiffness();


    print_vector(A);
    print_vector(B);
    print_vector(D);
    print_vector(C);


    return EXIT_SUCCESS;
}