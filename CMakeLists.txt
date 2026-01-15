cmake_minimum_required(VERSION 3.12)
project(global_polarization LANGUAGES CXX)

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------
option(BUILD_EXECUTABLES "Build standalone executables if sources have main()" OFF)
option(BUILD_BURN        "Build 'burn' runner: ./burn macro.C" ON)

# -----------------------------------------------------------------------------
# C++ settings
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -----------------------------------------------------------------------------
# ROOT
# -----------------------------------------------------------------------------
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree Physics MathCore GenVector Gpad)
message(STATUS "Using ROOT ${ROOT_VERSION} at ${ROOT_DIR}")

# ROOT CMake macros wrapper (your original pattern, but keep it early)
include(${ROOT_DIR}/modules/RootNewMacros.cmake OPTIONAL RESULT_VARIABLE _have_newmacros)
if(NOT _have_newmacros)
  include(${ROOT_DIR}/RootMacros.cmake)
endif()

# -----------------------------------------------------------------------------
# UniGen / MPDROOT paths
#   You had: set(UNIGEN_INCLUDE "{$MPDROOT}/include")  (wrong CMake syntax) :contentReference[oaicite:1]{index=1}
#   Here: allow either explicit -D..., or environment $MPDROOT.
# -----------------------------------------------------------------------------
set(UNIGEN_INCLUDE "" CACHE PATH "Path containing UParticle.h (UniGen headers)")
set(UNIGEN_LIBDIR  "" CACHE PATH "Path containing libUniGenFormat.so")

# convenient fallback to env MPDROOT if user didn't set cache vars
if((NOT UNIGEN_INCLUDE OR NOT UNIGEN_LIBDIR) AND DEFINED ENV{MPDROOT})
  set(_mpdroot "$ENV{MPDROOT}")
  if(NOT UNIGEN_INCLUDE)
    set(UNIGEN_INCLUDE "${_mpdroot}/include" CACHE PATH "" FORCE)
  endif()
  if(NOT UNIGEN_LIBDIR)
    set(UNIGEN_LIBDIR  "${_mpdroot}/lib" CACHE PATH "" FORCE)
  endif()
endif()

if(NOT UNIGEN_INCLUDE)
  message(FATAL_ERROR "UNIGEN_INCLUDE is not set. Provide -DUNIGEN_INCLUDE=/path/to/include (must contain UParticle.h)")
endif()
if(NOT UNIGEN_LIBDIR)
  message(FATAL_ERROR "UNIGEN_LIBDIR is not set. Provide -DUNIGEN_LIBDIR=/path/to/lib (must contain libUniGenFormat.so)")
endif()

message(STATUS "UniGen headers: ${UNIGEN_INCLUDE}")
message(STATUS "UniGen libdir : ${UNIGEN_LIBDIR}")

find_library(UNIGEN_LIB NAMES UniGenFormat PATHS ${UNIGEN_LIBDIR} REQUIRED)
message(STATUS "UniGenFormat found at: ${UNIGEN_LIB}")

# -----------------------------------------------------------------------------
# Sources from this repository (kept from your file)
# -----------------------------------------------------------------------------
set(GP_SOURCES
  src/calc_global_polarization.cpp
  src/calc_pol_vs_Nenh.cpp
  src/getters.cpp
  src/read_unigen_root.cpp
  src/set_lambda_parameterization.cpp
  src/simulate_lambda_decays.cpp
  src/add_enhanced_lambda.cpp
)

set(GP_INCLUDES
  include/calc_global_polarization.hpp
  include/calc_pol_vs_Nenh.hpp
  include/getters.hpp
  include/headers.hpp
  include/read_unigen_root.hpp
  include/set_lambda_parameterization.hpp
  include/simulate_lambda_decays.hpp
  include/add_enhanced_lambda.hpp
)

# -----------------------------------------------------------------------------
# Dictionary helper header + LinkDef (fixed + modern)
#   Your previous GPDictTypes.hh had a broken line: '#    include "UParticle.h"' :contentReference[oaicite:2]{index=2}
#   Also you duplicated __MAKECINT__/__CLING__ blocks. Use __ROOTCLING__ only.
# -----------------------------------------------------------------------------

# Dictionary files stored in repo
set(DICT_TYPES   ${CMAKE_CURRENT_SOURCE_DIR}/dict/GPDictTypes.hh)
set(DICT_LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/dict/GPDict_LinkDef.h)

# Extra include for rootcling to find UParticle.h
set(DICT_OPTS -I${UNIGEN_INCLUDE})

ROOT_GENERATE_DICTIONARY(G__GlobalPolarization
  ${DICT_TYPES}
  LINKDEF ${DICT_LINKDEF}
  OPTIONS ${DICT_OPTS}
)


# -----------------------------------------------------------------------------
# Build ONE shared library that contains:
#   - all your sources
#   - the generated dictionary source
#
# This replaces your "gp_dict MODULE + gp_macros MODULE + object trick" tangle :contentReference[oaicite:3]{index=3}
# -----------------------------------------------------------------------------
add_library(gp_macros SHARED
  ${GP_SOURCES}
  ${GP_INCLUDES}
  ${CMAKE_CURRENT_BINARY_DIR}/G__GlobalPolarization.cxx
)

target_include_directories(gp_macros PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${UNIGEN_INCLUDE}
)

target_link_libraries(gp_macros PUBLIC
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
  ROOT::Physics ROOT::MathCore ROOT::GenVector ROOT::Gpad
  ${UNIGEN_LIB}
)

# Optional: also build standalone executables if some sources contain main()
if(BUILD_EXECUTABLES)
  foreach(src ${GP_SOURCES})
    get_filename_component(name ${src} NAME_WE)
    add_executable(${name} ${src})
    target_include_directories(${name} PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      ${UNIGEN_INCLUDE}
    )
    target_link_libraries(${name} PRIVATE gp_macros)
  endforeach()
endif()

# -----------------------------------------------------------------------------
# 'burn' executable: ./burn macro.C
#   We generate a tiny runner in the build tree so you don't need to add a file.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 'burn' executable: ./burn macro.C
# -----------------------------------------------------------------------------
if(BUILD_BURN)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/burn.cpp
"#include <TROOT.h>\n"
"#include <stdexcept>\n"
"#include <iostream>\n"
"#include <string>\n"
"\n"
"int main(int n_args, char** args){\n"
"  if( n_args < 2 )\n"
"    throw std::runtime_error( \"No argumets provided\" );\n"
"  std::string macro{args[1]};\n"
"  std::string macro_full = macro+\"(\";\n"
"  for( int i=2; i<n_args; ++i ){\n"
"    macro_full += \"\\\"\" + std::string{args[i]} + \"\\\",\";\n"
"  }\n"
"  macro_full += ')';\n"
"  std::cout << \"Executing \" + macro_full << std::endl;\n"
"  gROOT->Macro(macro_full.c_str());\n"
"  return 0;\n"
"}\n"
  )

  add_executable(burn ${CMAKE_CURRENT_BINARY_DIR}/burn.cpp)
  target_link_libraries(burn PUBLIC gp_macros ROOT::Core ROOT::RIO)
endif()


# -----------------------------------------------------------------------------
# Optional: rootlogon.C (updated to load gp_macros only)
# -----------------------------------------------------------------------------
# file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/rootlogon.C
# "void rootlogon(){\n"
# "  gROOT->ProcessLine(\"gSystem->Load(\\\"${CMAKE_CURRENT_BINARY_DIR}/libgp_macros\\\")\");\n"
# "  gROOT->ProcessLine(\".include \\\"${CMAKE_CURRENT_SOURCE_DIR}/include\\\"\");\n"
# "}\n"
# )

# -----------------------------------------------------------------------------
# Install (optional)
# -----------------------------------------------------------------------------
include(GNUInstallDirs)
install(TARGETS gp_macros
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(BUILD_BURN)
  install(TARGETS burn
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
