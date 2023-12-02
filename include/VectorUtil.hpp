#ifndef VectorUtil_HPP
#define VectorUtil_HPP

#include <algorithm>
#include <fstream>
#include <vector>

#include <Eigen/Dense>

template <class T>
void write_vector(std::vector<T> vector, std::string filename)
{
    std::ofstream file(filename);

    if(file.is_open())
    {

    }

    std::cout << '{';

    for (std::size_t i = 0; i < vector.size() - 1; ++i)
        std::cout << vector[i] << ", ";

    std::cout << vector.back();

    std::cout << '}';
}

template <class T>
void print_vector(std::vector<T> vector)
{
    std::cout << '{';

    for (std::size_t i = 0; i < vector.size() - 1; ++i)
        std::cout << vector[i] << ", ";

    std::cout << vector.back();

    std::cout << '}';
}

template <class T>
inline void append_vector(std::vector<T> a, const std::vector<T> b)
{
    a.insert(std::end(a), std::begin(b), std::end(b));
}

// Check if vector a contains element x
template <class T>
inline bool vector_has_element(const std::vector<T> a, const T x)
{
    if (std::find(a.begin(), a.end(), x) != a.end())
        return true;

    else
        return false;
}

// Write Eigen matrix to CSV format file

// std::ofstream file("stiffness.txt");
// const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
//                                        Eigen::DontAlignCols,
//                                        ", ",
//                                        "\n");
// if (file.is_open())
//     file << element_stiffness.format(CSVFormat);
// file.close();

#endif // VectorUtil_HPP