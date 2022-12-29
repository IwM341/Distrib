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




add_executable(distrib_vector_grid distrib_vector_grid.cpp)
add_executable(matrix_mult matrix_modelling.cpp)



find_package(Boost COMPONENTS filesystem REQUIRED)
message("boost found:" ${Boost_FOUND})
message("BOOST_ROOT: " ${BOOST_ROOT})
message("boost include dir:" ${Boost_INCLUDE_DIR})
message("boost library dir:" ${Boost_LIBRARY_DIRS})
include_directories(distrib_vector_grid ${Boost_INCLUDE_DIR})
include_directories(matrix_mult ${Boost_INCLUDE_DIR})


find_package(OpenBLAS REQUIRED)

include_directories(matrix_mult ${OpenBLAS_INCLUDE_DIRS})


target_link_libraries(matrix_mult ${OpenBLAS_LIBRARIES})

option(TEST_CMD_PARSER "test parser" ON)
if(TEST_CMD_PARSER)
    add_executable(test_cmd_parser test_cmd_parser.cpp)
    include_directories(test_cmd_parser ${Boost_INCLUDE_DIR})
    target_link_libraries(test_cmd_parser  ${Boost_FILESYSTEM_LIBRARY})
endif(TEST_CMD_PARSER)
