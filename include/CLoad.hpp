#ifndef CLOAD_HPP
#define CLOAD_HPP

#include <vector>

class CLoad
{
public:
    CLoad(std::size_t node_tag,
          std::vector<bool> loading_dofs,
          std::vector<double> loading_values) : node_tag_(node_tag),
                                                loading_dofs_(loading_dofs),
                                                loading_values_(loading_values){};

    // The rule of five
    CLoad() = default;
    CLoad(CLoad const &) = default;
    CLoad &operator=(CLoad const &) = default;
    CLoad(CLoad &&) = default;
    CLoad &operator=(CLoad &&) = default;

    // Selectors
    inline std::size_t node_tag() const { return node_tag_; }
    std::vector<std::size_t> local_dofs() const;
    std::vector<std::size_t> global_dofs() const;
    std::vector<double> load_values() const;

private:
    std::size_t node_tag_{};
    std::vector<bool> loading_dofs_{};
    std::vector<double> loading_values_{};
};

#endif // CLOAD_HPP
