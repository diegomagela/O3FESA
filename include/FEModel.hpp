#ifndef FE_MODEL_H
#define FE_MODEL_H

#include <memory>
#include <map>
#include <unordered_map>

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
    FEModel(const std::string filename);

    // The rule of five
    FEModel() = delete;
    FEModel(FEModel const &) = default;
    FEModel &operator=(FEModel const &) = default;
    FEModel(FEModel &&) = default;
    FEModel &operator=(FEModel &&) = default;

    void print_nodes();
    void print_elements();
    void print_element_sets();
    void print_boundary();
    void print_cload();
    void print_dload();

private:
    /*
        Smart pointers for defining containers of derived classes and shared
        objects
    */

    typedef std::shared_ptr<Node> NodePtr;
    typedef std::unique_ptr<Element> ElementPtr;
    typedef std::shared_ptr<Boundary> BoundaryPtr;
    typedef std::shared_ptr<Material> MaterialPtr;
    typedef std::shared_ptr<Section> SectionPtr;

    /*
        TO TEST AND STUDY

        1) How to store all mesh data? Just store the data itself and spread
           this data for each specific class after reading it or read it and
           defined all classes at once (currently)?

        2) Use vector, map or unordered map?
    */

    // NODE TAG, NODE CLASS
    std::map<std::size_t, NodePtr> node_ptr_map{};

    // ELEMENT TAG, ELEMENT CLASSES
    std::map<std::size_t, ElementPtr> element_map{};

    /*
        Store then a list of element elements and access them through a map
        or store them as pointers (shared_ptr)?
    */

    // ELEMENT SET TAG, LIST OF ELEMENTS
    std::map<std::string, std::vector<std::size_t>> element_sets{};

    /*
        Set section and material in each ELEMENT object or
        get them from map when necessary?
    */

    // MATERIAL NAME, MATERIAL
    std::map<std::string, MaterialPtr> material_map{};

    // ELEMENT SET, SECTION
    std::map<std::string, SectionPtr> section_map{};

private:
    void add_node(const std::string &input);

    void add_element(const std::string &input,
                     const std::string &element_type,
                     std::vector<std::size_t> &element_set);

    void add_boundary(const std::string &input);

    void add_cload(const std::string &input);

    void add_dload(const std::string &input);

    void set_section_to_element();

};

#endif // FE_MODEL_H
