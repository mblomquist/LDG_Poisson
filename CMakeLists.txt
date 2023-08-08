cmake_minimum_required(VERSION 3.25)
project(LDG_Poisson)

set(CMAKE_CXX_STANDARD 17)

include_directories(/home/mblomquist/LAPACK/lapack-3.11.0/LAPACKE/include)

add_executable(LDG_Poisson src/main.cpp)
