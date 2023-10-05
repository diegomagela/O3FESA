# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Create a project
project(FEMSH)

# Set the project version number
project(FEMSH VERSION 0.1)

# Set the variable CMAKE_CXX_STANDARD to 20
# and the variable CMAKE_CXX_STANDARD_REQUIRED to True
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Use configure_file to configure and copy input file to binary directory
configure_file(${CMAKE_SOURCE_DIR}/mesh/mesh.inp.in mesh.inp)

# Add source files to the project
set(SOURCES
    src/Element.cpp
    src/IO.cpp
    src/Isotropic.cpp
    src/Lamina.cpp
    src/FEModelIO.cpp
    src/Section.cpp
    src/Q4.cpp
    src/Q9.cpp
    src/FEModel.cpp
    src/Node.cpp
    src/main.cpp
)

add_executable(femsh ${SOURCES})

# Use target_include_directories to include header files
target_include_directories(femsh PUBLIC "${CMAKE_SOURCE_DIR}/include")

# Include warning flags
# target_compile_options(femsh PRIVATE -Werror -Wall -Wextra)
target_compile_options(femsh PRIVATE -Wall -Wextra)


