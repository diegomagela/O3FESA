#ifndef FE_MODEL_IO_H
#define FE_MODEL_IO_H

#include <string>
#include <vector>

// Auxiliary IO functions to define the model from input file

bool find_keyword(const std::string &line, const std::string &keyword);

// Return element type
std::string get_element_type(const std::string &input);
// Return element set
std::string get_element_set(const std::string &input);
// Return material name
std::string get_material_name(const std::string &input);
// Return material type
std::string get_material_type(const std::string &input);
// Return section set name
std::string get_shell_section_set(const std::string &input);
// Return thermal expansion coefficient type
std::string get_expansion_type(const std::string &input); 
// Return step analysis type
std::string get_step_analysis_type(const std::string &input);



#endif // FE_MODEL_IO_H
