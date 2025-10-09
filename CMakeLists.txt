cmake_minimum_required(VERSION 3.23)
project(ShelfModel LANGUAGES Fortran C)

# ------------------------------------------------------------------------------
# Output structure
# ------------------------------------------------------------------------------
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)
endif()

# ------------------------------------------------------------------
# Per-config Fortran flags + optional heavy debug
# ------------------------------------------------------------------
# Let "debug" be fast-enough dev builds; "heavy-debug" adds extra checks.
option(HEAVY_DEBUG "Enable extra runtime checks/sanitizers in Debug" OFF)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Common
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -Wall -Wextra")

  # Debug (used by heavy-debug preset)
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan")

  # Normal dev mode (RelWithDebInfo): optimized + symbols
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")

  # Release: fast
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -funroll-loops -DNDEBUG")

  # MinSizeRel
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "-Os -DNDEBUG")

  # Extra: sanitizers in HEAVY_DEBUG (GNU & recent clang/gfortran)
  if(HEAVY_DEBUG)
    # Add to Debug only
    set(CMAKE_Fortran_FLAGS_DEBUG
        "${CMAKE_Fortran_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address,undefined,leak")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG
        "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address,undefined,leak")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # ifx/ifort
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -check all -traceback -fpe0 -ftrapuv")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -qopt-report=1 -DNDEBUG")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "-Os -DNDEBUG")
  # (Intel sanitizers differ; heavy-debug will still give strong runtime checks)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -N 1023")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -G2")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -G2 -DNDEBUG")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -h aggress -DNDEBUG")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "-Os -DNDEBUG")
endif()

# For multi-config generators (optional future-proofing)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/$<CONFIG>)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


# ------------------------------------------------------------------------------
# Optional libs
# ------------------------------------------------------------------------------
option(USE_OPENMP "Enable OpenMP" ON)

if(USE_OPENMP)
  find_package(OpenMP COMPONENTS Fortran)
endif()

# ------------------------------------------------------------------------------
# NetCDF (C + Fortran) detection 
# ------------------------------------------------------------------------------
option(USE_NETCDF "Enable NetCDF I/O support" ON)

if(USE_NETCDF)
  # Prefer CMake config packages provided by libnetcdf-dev & libnetcdff-dev
  find_package(NetCDF CONFIG)                 # C library; target: NetCDF::NetCDF
  find_package(NetCDF-Fortran CONFIG)         # Fortran;   target: NetCDF::NetCDF_Fortran

  if(TARGET NetCDF::NetCDF_Fortran)
    message(STATUS "Found NetCDF-Fortran (target NetCDF::NetCDF_Fortran)")
    # Link Fortran (which pulls C transitively on most distros)
    target_link_libraries(shelf_model PRIVATE NetCDF::NetCDF_Fortran)
  elseif(TARGET NetCDF::NetCDF)
    message(WARNING "Found NetCDF C but not Fortran. Install libnetcdff-dev or set USE_NETCDF=OFF.")
    target_link_libraries(shelf_model PRIVATE NetCDF::NetCDF)
  else()
    message(WARNING "NetCDF not found; continuing without it. Install libnetcdf-dev & libnetcdff-dev or set USE_NETCDF=OFF.")
  endIf()
endif()


# ------------------------------------------------------------------------------
# FABM submodule + embedded YAML
# ------------------------------------------------------------------------------

set(FABM_BASE "${CMAKE_SOURCE_DIR}/external/fabm")

if(NOT EXISTS "${FABM_BASE}/CMakeLists.txt")
  message(FATAL_ERROR
    "\nFABM submodule not found at: ${FABM_BASE}\n"
    "👉 Run:\n"
    "   git submodule update --init --recursive -- external/fabm\n"
    "Or set FABM_BASE to a valid path: -DFABM_BASE=/path/to/fabm\n"
  )
endif()

set(FABM_USE_YAML ON CACHE BOOL "" FORCE)
set(FABM_EMBED_YAFYAML ON CACHE BOOL "" FORCE)
set(FABM_DIMENSION 1 CACHE STRING "" FORCE)
set(FABM_VECTOR_DIM 1 CACHE STRING "" FORCE)
set(FABM_FORCED_HOST shelf_model)

# ---- Build FABM with warnings suppressed (only inside its subdir) ----
# Save current directory compile options (may be empty)
get_directory_property(_SAVED_COMPILE_OPTIONS COMPILE_OPTIONS)
if(NOT _SAVED_COMPILE_OPTIONS)
  set(_SAVED_COMPILE_OPTIONS "")
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set_directory_properties(PROPERTIES COMPILE_OPTIONS "-w")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set_directory_properties(PROPERTIES COMPILE_OPTIONS "-w")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  # adjust if needed; this often mutes most messages
  set_directory_properties(PROPERTIES COMPILE_OPTIONS "-hmsg=0")
endif()

# Use a single unique binary dir for FABM
add_subdirectory(${FABM_BASE} ${CMAKE_BINARY_DIR}/fabm-build)

# Restore options for the rest of your project (your code keeps -Wall -Wextra)
set_directory_properties(PROPERTIES COMPILE_OPTIONS "${_SAVED_COMPILE_OPTIONS}")



# ------------------------------------------------------------------------------
# Your model sources
# ------------------------------------------------------------------------------
file(GLOB_RECURSE MODEL_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_SOURCE_DIR}/src/*.f90"
  "${CMAKE_SOURCE_DIR}/src/*.F90"
)

add_executable(shelf_model ${MODEL_SOURCES})
set_property(TARGET shelf_model PROPERTY Fortran_STANDARD 2008)

# Pass FABM macros to host
target_compile_definitions(shelf_model PRIVATE
  FABM_NDIMS=${FABM_DIMENSION}
  FABM_VECTOR_DIM=${FABM_VECTOR_DIM}
)

# Include dirs for modules
target_include_directories(shelf_model PRIVATE
  ${CMAKE_Fortran_MODULE_DIRECTORY}
  ${CMAKE_SOURCE_DIR}/src
)

# Link libs
# (Swap 'fabm' → 'FABM::fabm' if your submodule exports namespaced targets)
target_link_libraries(shelf_model PRIVATE fabm)

if(TARGET OpenMP::OpenMP_Fortran)
  target_link_libraries(shelf_model PRIVATE OpenMP::OpenMP_Fortran)
endif()

if(USE_NETCDF AND NetCDF_FOUND)
  target_link_libraries(shelf_model PRIVATE NetCDF::NetCDF_Fortran)
endif()

message(STATUS "ShelfModel configured. FABM at: ${FABM_BASE}")

