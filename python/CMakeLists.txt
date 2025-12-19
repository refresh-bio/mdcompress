cmake_minimum_required(VERSION 3.15...3.27)
project(mdc)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()

find_package(Python 3.8 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/nanobind)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src/mdc_lib)


nanobind_add_module(mdc src/python_bind.cpp)

target_link_libraries(mdc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../bin/libmdc.a)

