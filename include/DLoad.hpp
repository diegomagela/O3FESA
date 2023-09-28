#ifndef DLOAD_HPP
#define DLOAD_HPP

#include <vector>

class DLoad
{
public:
    explicit DLoad(std::size_t element_tag,
                   bool loading,
                   double value) : element_tag_(element_tag),
                                   loading_(loading),
                                   value_(value){};

    // The rule of five
    DLoad() = default;
    DLoad(DLoad const &) = default;
    DLoad &operator=(DLoad const &) = default;
    DLoad(DLoad &&) = default;
    DLoad &operator=(DLoad &&) = default;

    bool has_load() { return loading_; };

private:
    std::size_t element_tag_{};
    bool loading_{};
    double value_{};
};

#endif // DLOAD_HPP
