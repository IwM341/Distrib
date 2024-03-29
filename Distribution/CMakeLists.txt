cmake_minimum_required(VERSION 3.5)

project(Distribution LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
#add_definitions(-DDELTAMASS)
#add_definitions(-DTEST_DISTRIBUTIONS)
#add_definitions(-DTEST_BORDERS)
add_definitions(-DTEST_SF_HISTO)
#add_definitions(-DTEST_TRAJECTORY)
#add_definitions(-DTEST_EQUIALENCE)
if(NOT DEFINED DM_PATH)
    set(DM_PATH "DMFramework\\utils")
endif()

if(NOT DEFINED GO_PATH)
    message("ERROR: set GO_PATH to GridObject headers")
endif()

include_directories(${DM_PATH})

include_directories(${GO_PATH})




#add_executable(distrib_vector_grid distrib_vector_grid_v1.cpp)
add_executable(distrib_vector_grid_grob distrib_vector_grid_grob.cpp)
add_executable(distrib_SD distrib_vector_grid_grob_SD.cpp)
#add_executable(intersec intersection.cpp)

add_executable(annihilation annihilation.cpp)
add_executable(annihilation_vector annihilation_vector.cpp)
add_executable(annihilation_rd annihilation_rd.cpp)
add_executable(annihilation_rd_vector annihilation_rd_vector.cpp)
add_executable(r_distrib r_distrib.cpp)

add_executable(matrix_mult matrix_modelling.cpp)
add_executable(matrix_mult_elastic matrix_modelling_elastic.cpp)
add_executable(matrix_mult_2 matrix2.cpp)

#add_executable(graphics graphs.cpp)
add_executable(conversions conversions.cpp)

#target_compile_definitions(graphics PUBLIC GNUPLOT_PATH="${GNUPLOT_PATH}")

message("required paths:")

message("BOOST_ROOT: " ${BOOST_ROOT})
message("OpenBLAS_DIR: " ${OpenBLAS_DIR})

find_package(Boost COMPONENTS filesystem)
message("boost found:" ${Boost_FOUND})
message("BOOST_ROOT: " ${BOOST_ROOT})
message("boost include dir:" ${Boost_INCLUDE_DIR})
message("boost library dir:" ${Boost_LIBRARY_DIRS})
message("boost filesystem lib:" ${Boost_FILESYSTEM_LIBRARY})

option(USE_OPENMP "use openmp" ON)
if(${USE_OPENMP})
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()
if(${Boost_FOUND})

#    include_directories(distrib_vector_grid ${Boost_INCLUDE_DIR})
#    include_directories(matrix_mult ${Boost_INCLUDE_DIR})
#    include_directories(graphics ${Boost_INCLUDE_DIR})

#    target_link_libraries(distrib_vector_grid PRIVATE ${Boost_FILESYSTEM_LIBRARY})
#    target_link_libraries(matrix_mult PRIVATE ${Boost_FILESYSTEM_LIBRARY})
#    target_link_libraries(graphics PRIVATE ${Boost_FILESYSTEM_LIBRARY})

 #   target_link_libraries(distrib_vector_grid PRIVATE Boost::filesystem)
    target_link_libraries(distrib_vector_grid_grob PRIVATE Boost::filesystem)
    target_link_libraries(distrib_SD PRIVATE Boost::filesystem)
    target_link_libraries(matrix_mult PRIVATE Boost::filesystem)
    #target_link_libraries(graphics PRIVATE Boost::filesystem)
    #target_link_libraries(intersec PRIVATE Boost::filesystem)
    target_link_libraries(conversions PRIVATE Boost::filesystem)
    target_link_libraries(annihilation PRIVATE Boost::filesystem)

    target_link_libraries(annihilation_vector PRIVATE Boost::filesystem)
    target_link_libraries(annihilation_rd PRIVATE Boost::filesystem)
    target_link_libraries(annihilation_rd_vector PRIVATE Boost::filesystem)
    target_link_libraries(r_distrib PRIVATE Boost::filesystem)
    target_link_libraries(matrix_mult_elastic PRIVATE Boost::filesystem)
    target_link_libraries(matrix_mult_2 PRIVATE Boost::filesystem)
endif()


find_package(OpenBLAS)
message("OpenBLAS found:" ${OpenBLAS_FOUND})
message("OpenBLAS_INCLUDE_DIRS: " ${OpenBLAS_INCLUDE_DIRS})
message("OpenBLAS_LIBRARIES: " ${OpenBLAS_LIBRARIES})


if(${OpenBLAS_FOUND})
    #include_directories(matrix_mult ${OpenBLAS_INCLUDE_DIRS})
    target_link_libraries(matrix_mult PRIVATE OpenBLAS::OpenBLAS)
    target_link_libraries(matrix_mult_elastic PRIVATE OpenBLAS::OpenBLAS)
    target_link_libraries(matrix_mult_2 PRIVATE OpenBLAS::OpenBLAS)
endif()


option(TEST_CMD_PARSER "test parser" OFF)
if(TEST_CMD_PARSER)
    add_executable(test_cmd_parser test_cmd_parser.cpp)
    if(${Boost_FOUND})
        target_link_libraries(test_cmd_parser PRIVATE Boost::filesystem)
        message("link to test_cmd_parser complete")
    endif()
endif(TEST_CMD_PARSER)

option(TEST_CMD_PARSER1 "test parser 1" OFF)
if(TEST_CMD_PARSER1)
    add_executable(test_cmd_parser1 test_cmd_parser1.cpp)
    if(${Boost_FOUND})
        target_link_libraries(test_cmd_parser1  PRIVATE Boost::filesystem)
        message("link to test_cmd_parser1 complete")
    endif()
endif(TEST_CMD_PARSER1)

option(TEST_SAVE_LOAD "test save load" OFF)
if(TEST_SAVE_LOAD)
    add_executable(test_save_load test_save_load.cpp)
    if(${Boost_FOUND})
        target_link_libraries(test_cmd_parser1  PRIVATE Boost::filesystem)
        message("link to test_save_load complete")
    endif()
endif(TEST_SAVE_LOAD)


option(PLOT_BM_FUNCTIONS "test save load" OFF)
if(PLOT_BM_FUNCTIONS)
    add_executable(plot_data plot_data.cpp)
    if(${Boost_FOUND})
        target_link_libraries(plot_data  PRIVATE Boost::filesystem)
        message("link to plot_data complete")
    endif()
endif()

option(PLOT_VELOCITY_DISTRIB "test save load" OFF)
if(PLOT_VELOCITY_DISTRIB)
    add_executable(plot_velocity plot_velocity.cpp)
    if(${Boost_FOUND})
        target_link_libraries(plot_velocity  PRIVATE Boost::filesystem)
        message("link to plot_data complete")
    endif()
endif()
