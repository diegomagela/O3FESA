#ifndef IO_HPP
#define IO_HPP

#include <string>
#include <vector>

// Auxiliary function to read the input file

// Remove string spaces
inline void remove_space(std::string &string);

// Split a string into a vector according to a specific separator
std::vector<std::string> split_string(const std::string input, const char sep);

// Print a splitted string (for debugging purposes)
void print_string_vector(const std::vector<std::string> string_vector);


#endif // IO_HPP
