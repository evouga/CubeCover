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
set(BUILD_VIZ_TOOLS OFF CACHE BOOL "Build Stuff for Volume Rendering")
set(BUILD_MINT_FRONTEND ON CACHE BOOL "Build the mint frame field solver front end.")
set(BUILD_MINT_UTILS ON CACHE BOOL "C++ utils used by some parts of mint, should be ported to mex")

set(USE_GUROBI ON CACHE BOOL "Use Gurobi for solving MIPs")
set(USE_COMISO OFF CACHE BOOL "Use Comiso for solving MIPs")

# Libigl
set(LIBIGL_WITH_COMISO $CACHE{USE_COMISO} CACHE BOOL "Build LIBIGL with CoMISo" FORCE)
include(libigl)

if(BUILD_SINGULARITY_VIEWER OR BUILD_ISOLINE_VIEWER OR BUILD_MINT_FRONTEND)
  # Find Polyscope
  find_package(Polyscope REQUIRED)
endif()

if(BUILD_VIZ_TOOLS)
  list(APPEND CMAKE_MODULE_PATH "/usr/local/lib/cmake/OpenVDB")
  find_package(OpenVDB REQUIRED)
endif()

if(USE_GUROBI)   
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
endif()

file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp)	  

add_library(${PROJECT_NAME} ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} igl::core)
if(USE_COMISO)
    target_link_libraries(${PROJECT_NAME} igl::comiso)
    target_compile_definitions(${PROJECT_NAME} PUBLIC HAS_COMISO)
endif()
if(USE_GUROBI)
    target_link_libraries(${PROJECT_NAME} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY} )
    target_compile_definitions(${PROJECT_NAME} PUBLIC HAS_GUROBI)
endif()

# Use C++14
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

install(TARGETS ${PROJECT_NAME})

# example programs

if(BUILD_SINGULARITY_VIEWER)
  add_subdirectory(tools/singularityviewer)
  add_subdirectory(tools/singularityviewer_fancy)
endif()

if(BUILD_ISOLINE_VIEWER)
  add_subdirectory(tools/isolineviewer)
  add_subdirectory(tools/isolineviewer_fancy)
  add_subdirectory(tools/boundarystreamlineviewer)
  add_subdirectory(tools/volumeparamviewer)
endif()

if(BUILD_VIZ_TOOLS)
  add_subdirectory(tools/hexextoOpenVDB)
endif()

if (BUILD_MINT_FRONTEND)
  add_subdirectory(tools/mint_frontend)
endif()

if (BUILD_MINT_UTILS)
  add_subdirectory(tools/proj_to_int_frames)
endif()



if(BUILD_BINARIES)
  add_subdirectory(tools/seamless)
  add_subdirectory(tools/integergrid)
endif()