#ifndef VectorUtil_HPP
#define VectorUtil_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace vct
{
    // Check if vector a contains element x
    template <class T>
    inline bool vector_has_element(const std::vector<T> &vector, const T x)
    {

        if (std::ranges::find(vector, x) != vector.end())
            return true;

        else
            return false;
    }

    // Extract elements of a vector a return a new vector with these elements
    template <class T>
    std::vector<T> extract_elements(const std::vector<T> &vector,
                                    const std::size_t first,
                                    const std::size_t last)
    {
        std::size_t size = last - first + 1; // number of elements
        std::vector<T> new_vector{};
        new_vector.reserve(size);

        for (std::size_t i = first; i <= last; ++i)
            new_vector.push_back(vector.at(i));

        return new_vector;
    }

    // Convert Eigen::VectorXd vector to std::vector<double>
    std::vector<double> eigen_vector_to_std(const Eigen::VectorXd &vector);

    // Convert a std::vector<double> to a Eigen::VectorXd
    inline Eigen::VectorXd
    std_to_eigen_vector(const std::vector<double> &vector)
    {
        return Eigen::Map<const Eigen::VectorXd>(vector.data(), vector.size());
    }

    // Write Eigen::_Xd to file
    void write_eigen(const Eigen::MatrixXd &object,
                     const std::string &filename);

} // namespace vecutl

#endif // VectorUtil_HPP