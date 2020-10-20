project("RPlus-Tree")

cmake_minimum_required(VERSION 3.15)

set(SFML_STATIC_LIBRARIES FALSE CACHE BOOL "Use static SFML librares")
if (SFML_STATIC_LIBRARIES)
	add_definitions(-DSFML_STATIC)
endif()
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/;${CMAKE_MODULE_PATH}")
find_package(SFML 2 COMPONENTS graphics window system)

if(SFML_FOUND)
	include_directories(${SFML_INCLUDE_DIR})
else()
	set(SFML_ROOT "" CACHE PATH "SFML top-level directory")
	message("\n-> SFML directory not found. Set SFML_ROOT to SFML's top-level path (containing \"include\" and \"lib\" directories).")
	message("-> Make sure the SFML libraries with the same configuration (Release/Debug, Static/Dynamic) exist.\n")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories($(CMAKE_CURRENT_SOURCE_DIR)/source)

target_link_libraries(rplus-tree_project ${SFML_LIBRARIES})

add_executable(rplus-tree_project source/main.cpp)