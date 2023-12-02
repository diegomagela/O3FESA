#ifndef TRIPLET_HPP
#define TRIPLET_HPP

#include <cstddef>
#include <iostream>

class Triplet
{
public:
    explicit Triplet(const std::size_t row,
                     const std::size_t col,
                     const double value) : row_(row),
                                           col_(col),
                                           value_(value){};

    // Selectors

    inline std::size_t row() const { return row_; }
    inline std::size_t col() const { return col_; }
    inline double value() const { return value_; }

    friend std::ostream &operator<<(std::ostream &os, const Triplet &triplet);

private:
    std::size_t row_;
    std::size_t col_;
    double value_;
};

#endif // TRIPLET_HPP
