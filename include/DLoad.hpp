#ifndef DLOAD_HPP
#define DLOAD_HPP

#include <vector>

class DLoad
{
public:
    explicit DLoad(std::size_t element_tag,
                   double value) : element_tag_(element_tag),
                                   value_(value){};

    // The rule of five
    DLoad() = default;
    DLoad(DLoad const &) = default;
    DLoad &operator=(DLoad const &) = default;
    DLoad(DLoad &&) = default;
    DLoad &operator=(DLoad &&) = default;

private:
    std::size_t element_tag_{};
    double value_{};
};

#endif // DLOAD_HPP
