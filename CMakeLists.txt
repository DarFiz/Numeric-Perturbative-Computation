cmake_minimum_required(VERSION 2.8.12.2)
project(NumericPerturbativeComputation)

# SET( CMAKE_CXX_COMPILER "/usr/bin/g++-5" )
# SET( CMAKE_CXX_COMPILER "/usr/bin/clang++-3.5" )
set(CMAKE_CXX_FLAGS "-g -O0 -Wall -Werror -std=c++11 -pthread")

#main
file(GLOB SOURCES "src/*.cpp")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
add_executable(main ${SOURCES})
target_include_directories(main PUBLIC ./include ./gnuplot/)
target_link_libraries (main armadillo)



