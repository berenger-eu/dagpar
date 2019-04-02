project(DAGPAR CXX)
cmake_minimum_required(VERSION 3.0)

# Set a default build type if none was specified
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Main -- Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#################################################
option(DAGPAR_USE_COVERAGE "Should be ON to enable code coverage" OFF)

if( DAGPAR_USE_COVERAGE )
    message(STATUS "Main -- Use code coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
    list(APPEND DAGPAR_LIBRARIES gcov)
    # Then one can use these commande:
    # lcov --directory CMakeFiles/ -c -o coverage-rapport.info
    # genhtml -o ./coverage -t "Coverage example" coverage-rapport.info
endif()

#################################################

if(METIS_DIR OR DEFINED ENV{METIS_DIR})
    message(STATUS "Main -- use metis")
    if(METIS_DIR)
        include_directories("${METIS_DIR}/include/")
        link_directories("${METIS_DIR}/lib/")
    else()
        include_directories("$ENV{METIS_DIR}/include/")
        link_directories("$ENV{METIS_DIR}/lib/")    
    endif()
    set(CMAKE_CXX_FLAGS "-DUSE_METIS")
    if(USE_SCOTCH_AS_METIS)
        if(SCOTCH_DIR)
            link_directories("${SCOTCH_DIR}/lib/")
        else()
            link_directories("$ENV{SCOTCH_DIR}/lib/")    
        endif() 
        list(APPEND DAGPAR_LIBRARIES scotchmetis scotch scotcherr pthread)
    endif()
    list(APPEND DAGPAR_LIBRARIES metis)
endif()

#################################################

option(DAGPAR_USE_ACYCLIC "Should be ON to enable code acyclic algo (need C++17)" OFF)
if(DAGPAR_USE_ACYCLIC)
    set(CMAKE_CXX_FLAGS "-DUSE_ACYCLIC")
endif()

#################################################

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")

SET(DAGPAR_EXTRA_FLAG "-m64 -march=native -Wextra -Wnon-virtual-dtor -Wshadow -Wpointer-arith -Wcast-qual -Wconversion  -Wall -Wno-sign-conversion -pedantic -Woverloaded-virtual -Wpointer-arith -Wcast-qual -Wconversion -Wno-error" CACHE STRING "Additional flags"  )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DAGPAR_EXTRA_FLAG}")

#################################################

message(STATUS "Main -- CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

file(GLOB_RECURSE cpp_for_lib ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
if(cpp_for_lib)
    file(GLOB_RECURSE hpp_for_lib ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)
    add_library(dagpar STATIC ${cpp_for_lib} ${hpp_for_lib})
endif()

#################################################


file(GLOB_RECURSE cpp_for_examples ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
foreach(an_example ${cpp_for_examples})
    get_filename_component(execname ${an_example} NAME_WE)
    
    message(STATUS "Main -- Add example ${execname} for ${an_example}")

    add_executable(${execname} ${an_example})
    if(cpp_for_lib)
        target_link_libraries(${execname} dagpar) 
    endif()
    if(DAGPAR_LIBRARIES)
        target_link_libraries(${execname} ${DAGPAR_LIBRARIES}) 
    endif()
endforeach()

#################################################
OPTION( DAGPAR_BUILD_UTESTS  "Set to ON to build the unit tests" ON )
if(DAGPAR_BUILD_UTESTS)
    enable_testing()

    file(GLOB_RECURSE cpp_for_utests ${CMAKE_CURRENT_SOURCE_DIR}/utests/*.cpp)
    foreach(an_example ${cpp_for_utests})
        get_filename_component(execname ${an_example} NAME_WE)
        
        message(STATUS "Main -- Add utests ${execname} for ${an_example}")

        add_executable(${execname} ${an_example})
        if(cpp_for_lib)
            target_link_libraries(${execname} dagpar) 
        endif()
        if(DAGPAR_LIBRARIES)
            target_link_libraries(${execname} ${DAGPAR_LIBRARIES}) 
        endif()
        add_test(${execname} ${CMAKE_CURRENT_BINARY_DIR}/${execname})
    endforeach()
endif()
