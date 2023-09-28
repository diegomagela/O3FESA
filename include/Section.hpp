#ifndef SECTION
#define SECTION

#include <vector>
#include <string>

class Section
{
public:
    Section(const std::vector<double> thickness,
            const std::vector<std::size_t> nip,
            const std::vector<std::string> material_name,
            const std::vector<int> orientation) : thickness_(thickness),
                                                  nip_(nip),
                                                  material_name_(material_name),
                                                  orientation_(orientation){};

    // The rule of five
    Section() = default;
    Section(Section const &) = default;
    Section &operator=(Section const &) = default;
    Section(Section &&) = default;
    Section &operator=(Section &&) = default;

private:
    std::vector<double> thickness_{};
    std::vector<std::size_t> nip_{}; // Number of integration points
    std::vector<std::string> material_name_{};
    std::vector<int> orientation_{};
};

#endif // SECTION
