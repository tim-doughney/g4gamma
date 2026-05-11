cmake_minimum_required(VERSION 3.18)
project(g4gamma LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Python is needed for pybind11
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# zlib for LARA tarball decompression (.tar.gz support)
find_package(ZLIB REQUIRED)

# pybind11 — prefer a system install if findable, else FetchContent
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    # GIT_TAG        v2.10.4
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
  )
  FetchContent_GetProperties(pybind11)

  if(NOT pybind11_POPULATED)
    message(STATUS "Fetching pybind11...")
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
  endif()
endif()

# Core library (pure C++, no Python dependency) so it can be reused
add_library(g4gamma_core STATIC
    src/DataPath.cc
    src/IsotopeKey.cc
    src/DecayData.cc
    src/PhotonEvap.cc
    src/FluorData.cc
    src/AugerData.cc
    src/EnsdfState.cc
    src/Geant4Provider.cc
    src/SandiaProvider.cc
    src/LaraProvider.cc
    src/ChainBuilder.cc
    src/Bateman.cc
    src/GammaSpectrum.cc
)
target_include_directories(g4gamma_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(g4gamma_core PRIVATE -Wall -Wextra -O2)
target_link_libraries(g4gamma_core PRIVATE ZLIB::ZLIB)

# pybind11 module
pybind11_add_module(g4gamma python/bindings.cc)
target_link_libraries(g4gamma PRIVATE g4gamma_core)
target_include_directories(g4gamma PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Optional: a tiny C++ CLI test driver
add_executable(g4gamma_cli src/cli.cc)
target_link_libraries(g4gamma_cli PRIVATE g4gamma_core)

install(TARGETS g4gamma DESTINATION .)

# ---------------------------------------------------------------------------
# SandiaDecay comparison test (requires third_party/SandiaDecay submodule)
# Run: git submodule update --init  then rebuild.
# ---------------------------------------------------------------------------
option(BUILD_SANDIA_COMPARE "Build sandia_compare Bateman validation binary" ON)
if(BUILD_SANDIA_COMPARE)
    set(SANDIA_DECAY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/SandiaDecay")
    if(EXISTS "${SANDIA_DECAY_DIR}/SandiaDecay.h")
        add_subdirectory(${SANDIA_DECAY_DIR})
        add_executable(sandia_compare test/sandia_compare.cc)
        target_link_libraries(sandia_compare PRIVATE g4gamma_core SandiaDecay)
        target_include_directories(sandia_compare PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${SANDIA_DECAY_DIR})
        message(STATUS "sandia_compare: enabled (SandiaDecay submodule found)")
    else()
        message(STATUS "sandia_compare: disabled (run: git submodule update --init)")
    endif()
endif()
