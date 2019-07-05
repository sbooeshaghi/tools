# cmake version to be used
cmake_minimum_required( VERSION 3.0 )

# project name
project( t2g )

# flags
SET(CMAKE_CXX_FLAGS "-std=c++11")
# files

# include
include_directories( ./include ./src )

# target
add_executable( t2g ./src/t2g_main.cpp )

# external libs

