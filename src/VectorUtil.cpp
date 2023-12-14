#ifndef VECTOR_UTIL_HPP
#define VECTOR_UTIL_HPP

#include "VectorUtil.hpp"

namespace vct
{
    std::vector<double> eigen_vector_to_std(const Eigen::VectorXd &vector)
    {
        std::vector<double> std_vector;
        std_vector.resize(vector.size());
        Eigen::VectorXd::Map(&std_vector[0], vector.size()) = vector;

        return std_vector;
    }

    void write_eigen(const Eigen::MatrixXd &object,
                     const std::string &filename)
    {
        std::ofstream file(filename);

        const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                               Eigen::DontAlignCols,
                                               ", ",
                                               "\n");
                                               
        if (file.is_open())
            file << object.format(CSVFormat);
        file.close();
    }
}

#endif // VECTOR_UTIL_HPP
