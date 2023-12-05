#ifndef VectorUtil_HPP
#define VectorUtil_HPP

#include <algorithm>
#include <fstream>
#include <vector>

#include <Eigen/Dense>

namespace vct
{
    template <class T>
    void write_vector(std::vector<T> vector, std::string filename)
    {
        std::ofstream file(filename);

        if (file.is_open())
        {

            std::cout << '{';

            for (std::size_t i = 0; i < vector.size() - 1; ++i)
                std::cout << vector[i] << ", ";

            std::cout << vector.back();

            std::cout << '}';
        }
    }

    template <class T>
    void print_vector(std::vector<T> vector)
    {
        std::cout << '{';

        for (std::size_t i = 0; i < vector.size() - 1; ++i)
            std::cout << vector[i] << ", ";

        std::cout << vector.back();

        std::cout << '}';

        std::cout << std::endl;
    }

    // Append vector B in the end of vector A
    template <class T>
    inline void append_vector(std::vector<T> &a, const std::vector<T> &b)
    {
        a.insert(a.end(), b.begin(), b.end());
    }

    // Check if vector a contains element x
    template <class T>
    inline bool vector_has_element(const std::vector<T> &a, const T x)
    {
        if (std::find(a.begin(), a.end(), x) != a.end())
            return true;

        else
            return false;
    }

    //
    template <class T>
    std::vector<T> extract_elements(const std::vector<T> &vector,
                                    const std::size_t first,
                                    const std::size_t last)
    {
        std::size_t size = last - first + 1; // number of elements
        std::vector<T> new_vector;
        new_vector.reserve(size);

        for (std::size_t i = first; i <= last; ++i)
            new_vector.push_back(vector.at(i));

        return new_vector;
    }
} // namespace vecutl

// FUNCTIONS TO CONSULT

// Write Eigen matrix to CSV format file

// std::ofstream file("stiffness.txt");
// const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
//                                        Eigen::DontAlignCols,
//                                        ", ",
//                                        "\n");
// if (file.is_open())
//     file << element_stiffness.format(CSVFormat);
// file.close();

// From std::vector to Eigen::VectorXd
// std::vector<double> A
// Eigen::VectorXd B
// B = Eigen::Map<Eigen::VectorXd>(A.data() , A.size());

#endif // VectorUtil_HPP