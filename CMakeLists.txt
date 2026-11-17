cmake_minimum_required(VERSION 3.16)
project(heat LANGUAGES CXX)

find_package(MPI COMPONENTS CXX)

add_executable(heat src/heat.cpp)
target_compile_features(heat PRIVATE cxx_std_17)

if(MPI_CXX_FOUND)
    target_link_libraries(heat PRIVATE MPI::MPI_CXX)
    target_compile_definitions(heat PRIVATE HEAT_HAVE_MPI=1)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(heat PRIVATE -O3 -ffp-contract=off -Wall -Wextra -Wpedantic)
endif()
