cmake_minimum_required(VERSION 2.8)
project(global_polarization LANGUAGES CXX)

# --- Options ---------------------------------------------------------------
option(BUILD_EXECUTABLES "Build standalone executables if sources have main()" OFF)

# --- C++ + warnings --------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- ROOT ------------------------------------------------------------------
# You need a ROOT 6 install on PATH/CMAKE_PREFIX_PATH.
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree Physics MathCore GenVector)
message(STATUS "Using ROOT ${ROOT_VERSION} at ${ROOT_DIR}")

# --- Where is UParticle.h? -------------------------------------------------
# Point this at your UniGen (MPDROOT) include folder that has UParticle.h.
# e.g. -DUNIGEN_INCLUDE=/path/to/mpdroot/simulation/generators/unigenFormat
# /mpdroot/simulation/generators/unigenFormat/UParticle.h
set(UNIGEN_INCLUDE "{$MPDROOT}/include")
if(UNIGEN_INCLUDE)
  message(STATUS "UniGen headers: ${UNIGEN_INCLUDE}")
else()
  message(WARNING "UNIGEN_INCLUDE is not set. Dictionary generation will fail if UParticle.h is not found.")
endif()
set(UNIGEN_LIBS "{$MPDROOT}/lib")

# --- Sources from this repository -----------------------------------------
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

set(GP_INCLUDES_DIR ${CMAKE_CURRENT_BINARY_DIR}/../include)

# --- Auto-generate a tiny header + LinkDef for the container dictionaries --
# We mimic gInterpreter->GenerateDictionary(...) calls.
set(DICT_DIR ${CMAKE_CURRENT_BINARY_DIR}/dict)
file(MAKE_DIRECTORY ${DICT_DIR})

file(WRITE ${DICT_DIR}/GPDictTypes.hh
"#pragma once\n"
"#include <vector>\n"
"#include \"TVector3.h\"\n"
"#include \"Math/Vector3D.h\"\n"
"#include \"UParticle.h\"\n"
"#ifdef __has_include\n"
"#  if __has_include(\"UParticle.h\")\n"
"#    include \"UParticle.h\"\n"
"#  endif\n"
"#endif\n"
)

file(WRITE ${DICT_DIR}/GPDict_LinkDef.h
"#ifdef __CLING__\n"
"#pragma clang diagnostic push\n"
"#pragma clang diagnostic ignored \"-Wkeyword-compat\"\n"
"#endif\n\n"
"#ifdef __MAKECINT__\n"
"#pragma link C++ class std::vector<UParticle>+;\n"
"#pragma link C++ class std::vector<TVector3>+;\n"
"#pragma link C++ class std::vector<ROOT::Math::XYZVector>+;\n"
"#endif\n\n"
"#ifdef __CLING__\n"
"#pragma link C++ class std::vector<UParticle>+;\n"
"#pragma link C++ class std::vector<TVector3>+;\n"
"#pragma link C++ class std::vector<ROOT::Math::XYZVector>+;\n"
"#pragma clang diagnostic pop\n"
"#endif\n"
)

# --- ROOT dictionary generation (modern CMake wrapper over rootcling) ------
include(${ROOT_DIR}/modules/RootNewMacros.cmake OPTIONAL RESULT_VARIABLE _have_newmacros)
if(NOT _have_newmacros)
  include(${ROOT_DIR}/RootMacros.cmake)  # fallback name on some installs
endif()

# Extra -I for UniGen headers if provided
set(DICT_OPTS "")
if(UNIGEN_INCLUDE)
  list(APPEND DICT_OPTS -I${UNIGEN_INCLUDE})
endif()

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/../rootlogon.C
"void rootlogon(){\n"
"    gROOT->ProcessLine(\"gSystem->Load(\\\"${CMAKE_CURRENT_BINARY_DIR}/libgp_dict\\\")\");\n"
"    gROOT->ProcessLine(\"gSystem->Load(\\\"${CMAKE_CURRENT_BINARY_DIR}/libgp_macros\\\")\");\n"
"    gROOT->ProcessLine(\".include \\\"${GP_INCLUDES_DIR}\\\"\");\n"
"    gROOT->ProcessLine(\"#include \\\"${GP_INCLUDES_DIR}/simulate_lambda_decays.hpp\\\"\");\n"
"}
"
)

# Build the dictionary as an object library we can link everywhere
add_library(gp_dict MODULE ${G__GlobalPolarization})
target_include_directories(gp_dict PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${UNIGEN_INCLUDE}
  /lhep/users/dflusova/mpdroot/install/include
)
#link_directories(/lhep/users/dflusova/mpdroot/install/lib)
find_library(UNIGEN_LIB NAMES UniGenFormat PATHS ${UNIGEN_LIBS} /lhep/users/dflusova/mpdroot/install/lib)

if(UNIGEN_LIB)
    message(STATUS "UniGenFormat found at: ${UNIGEN_LIB}")
else()
    message(FATAL_ERROR "UniGenFormat not found!")
endif()
target_link_libraries(gp_dict PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Physics ROOT::MathCore ROOT::GenVector ROOT::Gpad ${UNIGEN_LIB})

# Produce G__GlobalPolarization.cxx (plus pcm/rootmap in the build tree)
ROOT_GENERATE_DICTIONARY(G__GlobalPolarization
  MODULE gp_dict
  ${DICT_DIR}/GPDictTypes.hh
  LINKDEF ${DICT_DIR}/GPDict_LinkDef.h
  OPTIONS ${DICT_OPTS}
)

# For function exposure
# ROOT_GENERATE_DICTIONARY(G__GPFunctions
#   MODULE gp_funcs
#   ${DICT_DIR}/GPFunctions.hh  # Create empty header or include function headers
#   LINKDEF include/GPFunctions_LinkDef.h
#   OPTIONS ${DICT_OPTS}
# )

# --- Build a single MODULE (plugin) with your macros -----------------------
add_library(gp_macros MODULE ${GP_SOURCES} $<TARGET_OBJECTS:gp_dict> ${GP_INCLUDES})
target_include_directories(gp_macros PUBLIC ${UNIGEN_INCLUDE} /lhep/users/dflusova/mpdroot/install/include ${GP_INCLUDES_DIR})
target_link_libraries(gp_macros PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Physics ROOT::MathCore ROOT::GenVector ROOT::Gpad ${UNIGEN_LIB})

# Optional: also build standalone executables if your .cpp files contain main()
if(BUILD_EXECUTABLES)
  foreach(src ${GP_SOURCES})
    get_filename_component(name ${src} NAME_WE)
    add_executable(${name} ${src} $<TARGET_OBJECTS:gp_dict>)
    target_include_directories(${name} PRIVATE ${UNIGEN_INCLUDE} /lhep/users/dflusova/mpdroot/install/include ${GP_INCLUDES_DIR})
    target_link_libraries(${name} PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Physics ROOT::MathCore ROOT::GenVector ${UNIGEN_LIB})
  endforeach()
endif()

# add_executable(main main.cpp)
# target_include_directories(main PUBLIC include)
# # target_link_libraries(main PRIVATE some_library)
# target_link_libraries(main PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Physics ROOT::MathCore ROOT::GenVector ${UNIGEN_LIB})


# --- Install (optional) ----------------------------------------------------
include(GNUInstallDirs)
install(TARGETS gp_macros
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
