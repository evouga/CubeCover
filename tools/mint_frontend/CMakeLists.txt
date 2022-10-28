cmake_minimum_required(VERSION 3.11)
project(mint_frontend)

# Find Eigen
find_package(Eigen3 3.3 REQUIRED)

file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp)	  

add_executable(${PROJECT_NAME} ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR} ../../src)

include(FetchContent)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
FetchContent_MakeAvailable(json)

# Use C++14
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)


IF (WIN32)
    target_link_libraries(${PROJECT_NAME} cubecover polyscope nlohmann_json::nlohmann_jsons)
ELSE()
    target_link_libraries(${PROJECT_NAME} cubecover polyscope nlohmann_json::nlohmann_json stdc++fs)
ENDIF()


# stdc++fs
