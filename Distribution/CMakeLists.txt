cmake_minimum_required(VERSION 3.5)

project(Distribution LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#add_definitions(-DDELTAMASS)
#add_definitions(-DTEST_DISTRIBUTIONS)
#add_definitions(-DTEST_BORDERS)
add_definitions(-DTEST_SF_HISTO)
#add_definitions(-DTEST_TRAJECTORY)
#add_definitions(-DTEST_EQUIALENCE)
include_directories("..\\..\\DMFramework\\utils" )

add_executable(Distribution main_2Comp.cpp)
