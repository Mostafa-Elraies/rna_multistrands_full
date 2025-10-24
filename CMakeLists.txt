cmake_minimum_required(VERSION 3.10)
project(MyViennaRNAProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Ensure Conda environment variables are set
if(NOT DEFINED ENV{CONDA_PREFIX})
    message(FATAL_ERROR "Conda environment is not activated. Please activate your environment.")
endif()

# Set ViennaRNA paths from Conda
set(VIENNARNA_INCLUDE_DIR $ENV{CONDA_PREFIX}/include)
set(VIENNARNA_LIBRARY $ENV{CONDA_PREFIX}/lib/libRNA.a)

# Include directories
include_directories(${VIENNARNA_INCLUDE_DIR})

# Add the executable
add_executable(main nussinov.cpp global_variables.cpp utilities.cpp)

# Define the NUSSINOV_TEST macro
target_compile_definitions(main PRIVATE NUSSINOV_TEST)

# Link against the ViennaRNA library
target_link_libraries(main PRIVATE ${VIENNARNA_LIBRARY} m)

