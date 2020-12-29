cmake_minimum_required(VERSION 3.1)
project(cubecover)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(MSVC)
  set(MSVC_RUNTIME "dynamic")
  set(BLA_STATIC ON)
endif()

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

# libigl
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          OFF)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  OFF)
option(LIBIGL_WITH_XML               "Use XML"            OFF)

find_package(LIBIGL REQUIRED)
message("found libigl at ${LIBIGL_INCLUDE_DIR}")

# use gurobi
find_package(GUROBI)
if(NOT GUROBI_FOUND)
  message("Unable to find Gurobi. Make sure that the env GUROBI_HOME is set")
  set (GUROBI_INCLUDE_DIRS " C:/gurobi901/win64/include")
  set (GUROBI_CXX_LIBRARY "C:/gurobi901/win64/lib/gurobi_c++md2019.lib")
  set (GUROBI_LIBRARY "C:/gurobi901/win64/lib/gurobi90.lib")
else()
  message("Found Gurobi in : " ${GUROBI_INCLUDE_DIRS})
  message("Found Gurobi cxx lib in : " ${GUROBI_CXX_LIBRARY})
  message("Found Gurobi lib in : " ${GUROBI_LIBRARY})
  include_directories(${GUROBI_INCLUDE_DIRS})
endif()

# Add your project files
file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp)	  

add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})

# Use C++14
set_target_properties(${PROJECT_NAME}_bin PROPERTIES CXX_STANDARD 14)
set_target_properties(${PROJECT_NAME}_bin PROPERTIES CXX_STANDARD_REQUIRED ON)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME}_bin OpenMP::OpenMP_CXX)
endif()