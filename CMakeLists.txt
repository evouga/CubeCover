cmake_minimum_required(VERSION 3.1)
project(cubecover)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(BUILD_SINGULARITY_VIEWER ON CACHE BOOL "Build Singularity Viewer")
set(BUILD_ISOLINE_VIEWER ON CACHE BOOL "Build Isoline Viewer")
set(BUILD_BINARIES ON CACHE BOOL "Build Command-line Binaries")

# Find Eigen
find_package(Eigen3 3.3 REQUIRED)

if(BUILD_SINGULARITY_VIEWER OR BUILD_ISOLINE_VIEWER)
  # Find Polyscope
  find_package(Polyscope REQUIRED)
endif()

# use gurobi
find_package(GUROBI)
if(NOT GUROBI_FOUND)
  message(SEND_ERROR "Unable to find Gurobi. Make sure that the env GUROBI_HOME is set")
else()
  message("Found Gurobi in : " ${GUROBI_INCLUDE_DIRS})
  message("Found Gurobi cxx lib in : " ${GUROBI_CXX_LIBRARY})
  message("Found Gurobi lib in : " ${GUROBI_LIBRARY})
  include_directories(${GUROBI_INCLUDE_DIRS})
endif()

file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp)	  

add_library(${PROJECT_NAME} ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})

# Use C++14
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

install(TARGETS ${PROJECT_NAME})

# example programs

if(BUILD_SINGULARITY_VIEWER)
  add_subdirectory(tools/singularityviewer)
endif()

if(BUILD_ISOLINE_VIEWER)
  add_subdirectory(tools/isolineviewer)
endif()

if(BUILD_BINARIES)
  add_subdirectory(tools/seamless)
  add_subdirectory(tools/integergrid)
endif()