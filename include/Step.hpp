#ifndef STEP_HPP
#define STEP_HPP

#include <string>
#include <map>
#include <memory>
#include <vector>

#include "Boundary.hpp"
#include "CLoad.hpp"
#include "DLoad.hpp"

class Step
{
public:
    typedef std::shared_ptr<Boundary> BoundaryPtr;
    typedef std::shared_ptr<CLoad> CLoadPtr;
    typedef std::shared_ptr<DLoad> DLoadPtr;

    Step(const std::string &type,
         const std::map<std::size_t, BoundaryPtr> boundary_map,
         const std::map<std::size_t, CLoadPtr> cload_map,
         const std::map<std::size_t, DLoadPtr> dload_map) : type_(type),
                                                            boundary_map_(boundary_map),
                                                            cload_map_(cload_map),
                                                            dload_map_(dload_map){};
    ~Step(){};

private:
    const std::string type_{};
    const std::map<std::size_t, BoundaryPtr> boundary_map_{};
    const std::map<std::size_t, CLoadPtr> cload_map_{};
    const std::map<std::size_t, DLoadPtr> dload_map_{};
};

#endif // STEP_HPP
