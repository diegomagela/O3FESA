#include "IO.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>

std::vector<std::string> split_string(const std::string input, const char sep)
{
    std::stringstream ss(input);
    std::vector<std::string> splitted_string;

    while (ss.good())
    {
        std::string substr;
        std::getline(ss, substr, sep);
        // substr = remove_space(substr);
        remove_space(substr);
        splitted_string.push_back(substr);
    }

    return splitted_string;
};

void print_string_vector(const std::vector<std::string> string_vector)
{
    std::cout << "Vector size: " << string_vector.size() << '\n';

    for (const auto &string : string_vector)
        std::cout << string << '\n';
}
