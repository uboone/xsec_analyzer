# Require CMake 3.16 or higher
# 3.16 was released in 2019 and includes:
# - Precompiled header support (needed for heavy single-header dependencies)
# - Modern target propagation features
# - Unity builds support
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(XSecAnalyzer VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard (C++17 required for XSecAnalyzer and YODEL)
# Later standards (C++20, C++23) are fine and will be used if available
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add our include directory to the include path
include_directories(${CMAKE_SOURCE_DIR}/include)

# -----------------------------------------------------------------------------
# Find ROOT installation
# -----------------------------------------------------------------------------
# Option to force use of root-config even if find_package would work
option(FORCE_ROOT_CONFIG "Force use of root-config instead of find_package" OFF)

# Try the modern FindROOT.cmake first (available in ROOT 6.16+)
# If that fails, fall back to using root-config directly
if(NOT FORCE_ROOT_CONFIG)
  set(ROOT_REQUIRED_COMPONENTS Core RIO Tree Hist Graf Gpad Physics TreePlayer Rint)
  
  # Temporarily suppress warnings during ROOT detection
  # ROOT's CMake config warns about every optional library not installed
  # We'll manually report any missing components we actually need
  set(_cmake_message_log_level ${CMAKE_MESSAGE_LOG_LEVEL})
  set(CMAKE_MESSAGE_LOG_LEVEL ERROR)
  
  find_package(ROOT 6.0 CONFIG QUIET)
  
  # Restore normal message level
  set(CMAKE_MESSAGE_LOG_LEVEL ${_cmake_message_log_level})
  
  # Check if we got all the components we need
  if(ROOT_FOUND)
    set(ROOT_HAS_ALL_COMPONENTS TRUE)
    set(ROOT_MISSING_COMPONENTS "")
    
    foreach(comp ${ROOT_REQUIRED_COMPONENTS})
      if(NOT TARGET ROOT::${comp})
        set(ROOT_HAS_ALL_COMPONENTS FALSE)
        list(APPEND ROOT_MISSING_COMPONENTS ${comp})
      endif()
    endforeach()
    
    if(NOT ROOT_HAS_ALL_COMPONENTS)
      set(ROOT_FOUND FALSE)
      # Provide a clear, actionable message about what's actually missing
      message(WARNING "Required ROOT components not found: ${ROOT_MISSING_COMPONENTS}")
      message(STATUS "Will use root-config fallback instead")
    endif()
  endif()
endif()

if(NOT ROOT_FOUND OR FORCE_ROOT_CONFIG)
  # Fallback: use root-config executable directly
  # This is more fragile but works with older ROOT 6 versions
  find_program(ROOT_CONFIG_EXECUTABLE root-config)
  
  if(NOT ROOT_CONFIG_EXECUTABLE)
    message(FATAL_ERROR 
      "Could not find a working ROOT installation.\n"
      "Please ensure that the root-config executable is on your PATH and try again.")
  endif()
  
  # Get ROOT configuration using root-config
  execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --cflags 
                  OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --libs 
                  OUTPUT_VARIABLE ROOT_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --incdir
                  OUTPUT_VARIABLE ROOT_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --version
                  OUTPUT_VARIABLE ROOT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
  
  if(FORCE_ROOT_CONFIG)
    message(STATUS "Found ROOT ${ROOT_VERSION} via root-config (forced)")
  else()
    message(STATUS "Found ROOT ${ROOT_VERSION} via root-config (fallback - some CMake components unavailable)")
  endif()
  
  # Convert the flags to lists for CMake
  separate_arguments(ROOT_CXX_FLAGS)
  separate_arguments(ROOT_LIBRARIES)
  
  # Add ROOT include directories
  include_directories(${ROOT_INCLUDE_DIRS})
  
  # Find rootcling for dictionary generation
  find_program(ROOT_CLING_EXECUTABLE rootcling)
  if(NOT ROOT_CLING_EXECUTABLE)
    message(FATAL_ERROR "Could not find rootcling executable")
  endif()
else()
  message(STATUS "Found ROOT ${ROOT_VERSION} via CMake config")
endif()

# -----------------------------------------------------------------------------
# Build configuration
# -----------------------------------------------------------------------------
# Set optimization flags based on build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_compile_options(-O0 -g)
else()
  add_compile_options(-O3)
endif()

# -----------------------------------------------------------------------------
# Generate ROOT dictionary
# -----------------------------------------------------------------------------
set(DICT_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/dictionaries.cxx)

if(ROOT_FOUND AND TARGET ROOT::Core)
  # Modern ROOT with CMake support - use the provided function
  # Note: ROOT_GENERATE_DICTIONARY is provided by ROOTConfig.cmake
  include(${ROOT_USE_FILE})
  ROOT_GENERATE_DICTIONARY(G__dictionaries 
    LINKDEF ${CMAKE_SOURCE_DIR}/LinkDef.hh
    MODULE XSecAnalyzer
  )
  set(DICT_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/G__dictionaries.cxx)
else()
  # Fallback: call rootcling directly
  add_custom_command(
    OUTPUT ${DICT_OUTPUT}
    COMMAND ${ROOT_CLING_EXECUTABLE} -f ${DICT_OUTPUT} -c ${CMAKE_SOURCE_DIR}/LinkDef.hh
    DEPENDS ${CMAKE_SOURCE_DIR}/LinkDef.hh
    COMMENT "Generating ROOT dictionary"
  )
  set(DICT_SOURCE ${DICT_OUTPUT})
endif()

# -----------------------------------------------------------------------------
# Build the main shared library
# -----------------------------------------------------------------------------
# Collect all source files
# CONFIGURE_DEPENDS tells CMake to check for new/removed files and reconfigure if needed
file(GLOB_RECURSE SHARED_SOURCES CONFIGURE_DEPENDS
  src/binning/*.cxx
  src/selections/*.cxx
  src/utils/*.cxx
)

# Create the shared library
add_library(XSecAnalyzer SHARED
  ${SHARED_SOURCES}
  ${DICT_SOURCE}
)

# Link against ROOT libraries
if(ROOT_FOUND AND TARGET ROOT::Core)
  # Modern ROOT: use imported targets
  # Specify all the ROOT components we actually use to avoid linking issues
  # with incomplete ROOT installations (e.g., Fedora's modular packages)
  target_link_libraries(XSecAnalyzer PUBLIC 
    ROOT::Core 
    ROOT::RIO 
    ROOT::Tree 
    ROOT::Hist
    ROOT::Graf      # TLine
    ROOT::Gpad      # TCanvas, TPad, TLegend
    ROOT::Physics   # TVector3, TVector2, TLorentzVector
    ROOT::TreePlayer # TTreeFormula
    ROOT::Rint      # TRint (for binscheme app)
  )
else()
  # Fallback: use the libraries from root-config
  target_link_libraries(XSecAnalyzer PUBLIC ${ROOT_LIBRARIES})
  target_compile_options(XSecAnalyzer PUBLIC ${ROOT_CXX_FLAGS})
endif()

# Precompile heavy headers to speed up compilation
# This is especially important for single-header libraries like fkYAML
target_precompile_headers(XSecAnalyzer PUBLIC
  # Standard library headers used throughout
  <string>
  <vector>
  <memory>
  <map>
  # Add your single-header dependencies here when ready:
  # <fkYAML/node.hpp>  # ~15k lines, very heavy
  # <yodel/yodel.hpp>
)

# Set library output directory to build/lib
set_target_properties(XSecAnalyzer PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# -----------------------------------------------------------------------------
# Build executables
# -----------------------------------------------------------------------------
# Helper function to create executables that link against our library
function(add_xsec_executable target_name source_file)
  add_executable(${target_name} ${source_file})
  target_link_libraries(${target_name} PRIVATE XSecAnalyzer)
  
  # Set executable output directory to build/bin
  set_target_properties(${target_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
  )
  
  # Executables automatically inherit the precompiled headers from XSecAnalyzer
  
  # For NuMI apps, force O3 optimization even in debug builds
  if(source_file MATCHES "NuMI")
    target_compile_options(${target_name} PRIVATE -O3)
  endif()
endfunction()

# Build all executables
add_xsec_executable(ProcessNTuples src/app/ProcessNTuples.C)
add_xsec_executable(univmake src/app/univmake.C)
add_xsec_executable(SlicePlots src/app/Slice_Plots.C)
add_xsec_executable(Unfolder src/app/Unfolder.C)
add_xsec_executable(UnfolderNuMI src/app/UnfolderNuMI.C)
add_xsec_executable(BinScheme src/app/binscheme.C)
add_xsec_executable(StandaloneUnfold src/app/standalone_unfold.C)
add_xsec_executable(AddFakeWeights src/app/NuMI/addFakeWeights.cpp)
add_xsec_executable(AddBeamlineGeometryWeights src/app/NuMI/addBeamlineGeometryWeightsToMap.cpp)

# -----------------------------------------------------------------------------
# Handle script files (xsroot and xsnotebook)
# -----------------------------------------------------------------------------
# Copy shell scripts to build/bin directory
configure_file(src/app/xsroot ${CMAKE_BINARY_DIR}/bin/xsroot COPYONLY)
configure_file(src/app/xsnotebook ${CMAKE_BINARY_DIR}/bin/xsnotebook COPYONLY)

# Make scripts executable (CMake 3.19+ has better syntax, but this works for 3.16+)
file(CHMOD ${CMAKE_BINARY_DIR}/bin/xsroot 
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(CHMOD ${CMAKE_BINARY_DIR}/bin/xsnotebook
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# -----------------------------------------------------------------------------
# Installation rules (optional, but good practice)
# -----------------------------------------------------------------------------
# Note: Most users will use directly from the build directory via the setup script
install(TARGETS XSecAnalyzer
  LIBRARY DESTINATION lib
)

install(TARGETS 
  ProcessNTuples univmake SlicePlots Unfolder UnfolderNuMI 
  BinScheme StandaloneUnfold AddFakeWeights AddBeamlineGeometryWeights
  RUNTIME DESTINATION bin
)

install(PROGRAMS 
  ${CMAKE_BINARY_DIR}/bin/xsroot 
  ${CMAKE_BINARY_DIR}/bin/xsnotebook
  DESTINATION bin
)
