cmake_minimum_required(VERSION 3.14)
project(BoundaryConditions)

set(CMAKE_CXX_STANDARD 14)
include_directories(C:/dealii-9.0.1/dealii-9.0.1/include)


add_executable(BoundaryConditions main.cpp forces_inclusions.h membrane.cpp assemble.cpp solve.cpp setup.cpp output.cpp simulate_membrane.cpp)