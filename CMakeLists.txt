cmake_minimum_required(VERSION 3.16)

project(heat LANGUAGES CXX)

# Find MPI
find_package(MPI REQUIRED COMPONENTS CXX)

# Build the executable
add_executable(heat src/heat.cpp)

# C++17
target_compile_features(heat PRIVATE cxx_std_17)

# Connect MPI to the executable
target_link_libraries(heat PRIVATE MPI::MPI_CXX)

# Optimization and numerical reproducibility
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(heat PRIVATE
        -O3
        -ffp-contract=off
        -Wall
        -Wextra
        -Wpedantic
    )
endif()
