#ifndef FE_MODEL_H
#define FE_MODEL_H

#include <string>
#include <fstream>
#include <memory>
#include <map>
#include <map>

#include "Element.hpp"
#include "Material.hpp"
#include "Node.hpp"
#include "Section.hpp"

/*
    TODO

    1) For better performance, use maps to store nodes and elements.
       - I'm not sure if this will indeed result in better performance.

    2) For improved performance, use an element's tag when storing element sets.
       - Solved (?): Using pointers to assign nodes to elements.

    3) Currently, the input file must follow a specific order of keywords, but the code
       should be able to read them in any order.
       - Is it a good idea to reopen the file multiple times to read a specific keyword?
*/

// Class that owns all model properties
class FEModel
{
public:
    FEModel(const std::string filename) : filename_(filename)
    {
        std::fstream input(filename_);
        std::string line;

        if (!input.is_open())
        {
            std::cerr << "Input file not found!" << std::endl;
            abort();
        }
    };

    // The rule of five
    FEModel() = delete;
    FEModel(FEModel const &) = default;
    FEModel &operator=(FEModel const &) = default;
    FEModel(FEModel &&) = default;
    FEModel &operator=(FEModel &&) = default;

    // Read mesh
    void read_input();

    void print_nodes();
    void print_elements();

private:
    // Smart pointers for defining containers of derived classes and shared
    // objects

    typedef std::shared_ptr<Boundary> BoundaryPtr;
    typedef std::shared_ptr<CLoad> CLoadPtr;
    typedef std::shared_ptr<Node> NodePtr;
    typedef std::shared_ptr<Material> MaterialPtr;
    typedef std::shared_ptr<Section> SectionPtr;
    typedef std::shared_ptr<DLoad> DLoadPtr;
    typedef std::shared_ptr<Element> ElementPtr;

private:
    const std::string filename_{};

    /*
        TO TEST AND STUDY

        1) How to store all mesh data? Just store the data itself and spread
           this data for each specific class after reading it or read it and
           defined all classes at once (currently)?

        2) Use vector, map or unordered map? Need some benchmarks, but unordered
           map should be the best option here.
    */

    // NODE TAG, BOUNDARY
    std::map<std::size_t, BoundaryPtr> boundary_map{};

    // NODE TAG, CLOAD
    std::map<std::size_t, CLoadPtr> cload_map{};

    // NODE TAG, NODE CLASS
    std::map<std::size_t, NodePtr> node_map{};

    // MATERIAL NAME, MATERIAL CLASS
    std::map<std::string, MaterialPtr> material_map{};

    // ELEMENT SET, SECTION
    std::map<std::string, SectionPtr> section_map{};

    // ELEMENT TAG, DLOAD
    std::map<std::size_t, DLoadPtr> dload_map{};

    // ELEMENT TAG, ELEMENT CLASSES
    std::map<std::size_t, ElementPtr> element_map{};

    void read_boundary();
    void read_cload();
    void read_nodes();
    void read_materials();
    void read_sections();
    void read_dload();
    void read_elements();
};

#endif // FE_MODEL_H
