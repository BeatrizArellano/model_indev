# Tries CONFIG -> nf-config -> pkg-config -> ENV fallbacks
# Exposes imported targets and classic variables.
include_guard(GLOBAL)

# --- 1) Try CMake config packages (modern) ---
set(_nc_ok FALSE)
find_package(NetCDF-Fortran CONFIG QUIET)
if (TARGET NetCDF::NetCDF_Fortran)
  set(_nc_ok TRUE)
  # Ensure C target exists too (some configs do not export it)
  if (NOT TARGET NetCDF::NetCDF)
    find_package(NetCDF CONFIG QUIET)
  endif()
endif()

# --- 2) nf-config (Linux/macOS, most HPC) ---
if (NOT _nc_ok)
  find_program(NetCDF_CONFIG_EXECUTABLE NAMES nf-config DOC "NetCDF-Fortran config tool")
  if (NetCDF_CONFIG_EXECUTABLE)
    execute_process(COMMAND ${NetCDF_CONFIG_EXECUTABLE} --includedir
                    OUTPUT_VARIABLE _NC_F_INC OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${NetCDF_CONFIG_EXECUTABLE} --flibs
                    OUTPUT_VARIABLE _NC_F_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
    # Optional: also get C libs (some nf-config provide --clibs)
    execute_process(COMMAND ${NetCDF_CONFIG_EXECUTABLE} --clibs
                    OUTPUT_VARIABLE _NC_C_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    set(_nc_ok TRUE)

    # Create imported targets from nf-config info
    if (NOT TARGET NetCDF::NetCDF_Fortran)
      add_library(NetCDF::NetCDF_Fortran INTERFACE IMPORTED)
      target_include_directories(NetCDF::NetCDF_Fortran INTERFACE ${_NC_F_INC})
      separate_arguments(_NC_F_LIBS) # split "-L... -lnetcdff -lnetcdf -lhdf5 ..."
      target_link_libraries(NetCDF::NetCDF_Fortran INTERFACE ${_NC_F_LIBS})
    endif()
    if (NOT TARGET NetCDF::NetCDF AND _NC_C_LIBS)
      add_library(NetCDF::NetCDF INTERFACE IMPORTED)
      separate_arguments(_NC_C_LIBS)
      target_link_libraries(NetCDF::NetCDF INTERFACE ${_NC_C_LIBS})
    endif()

    # Classic variables for compatibility
    set(NetCDF_INCLUDE_DIRS ${_NC_F_INC})
    set(NetCDF_LIBRARIES    ${_NC_F_LIBS})
    set(NetCDF_Fortran_FOUND TRUE)
  endif()
endif()

# --- 3) pkg-config (Linux/macOS) ---
if (NOT _nc_ok)
  find_package(PkgConfig QUIET)
  if (PkgConfig_FOUND)
    pkg_check_modules(NETCDF_F QUIET IMPORTED_TARGET netcdf-fortran)
    pkg_check_modules(NETCDF_C QUIET IMPORTED_TARGET netcdf)

    if (NETCDF_F_FOUND)
      add_library(NetCDF::NetCDF_Fortran INTERFACE IMPORTED)
      target_include_directories(NetCDF::NetCDF_Fortran INTERFACE ${NETCDF_F_INCLUDE_DIRS} ${NETCDF_C_INCLUDE_DIRS})
      target_link_libraries(NetCDF::NetCDF_Fortran INTERFACE ${NETCDF_F_LINK_LIBRARIES} ${NETCDF_C_LINK_LIBRARIES})
      set(NetCDF_INCLUDE_DIRS ${NETCDF_F_INCLUDE_DIRS} ${NETCDF_C_INCLUDE_DIRS})
      set(NetCDF_LIBRARIES    ${NETCDF_F_LINK_LIBRARIES} ${NETCDF_C_LINK_LIBRARIES})
      set(_nc_ok TRUE)
      set(NetCDF_Fortran_FOUND TRUE)
    endif()

    if (NETCDF_C_FOUND AND NOT TARGET NetCDF::NetCDF)
      add_library(NetCDF::NetCDF INTERFACE IMPORTED)
      target_include_directories(NetCDF::NetCDF INTERFACE ${NETCDF_C_INCLUDE_DIRS})
      target_link_libraries(NetCDF::NetCDF INTERFACE ${NETCDF_C_LINK_LIBRARIES})
    endif()
  endif()
endif()

# --- 4) ENV fallbacks / manual hints (last resort, also helps Windows) ---
if (NOT _nc_ok)
  # User can provide NETCDFINC, NETCDFLIBDIR, NETCDFLIBNAME
  if (DEFINED ENV{NETCDFINC})
    set(NetCDF_INCLUDE_DIRS $ENV{NETCDFINC})
  endif()
  set(_LIBS "")
  if (DEFINED ENV{NETCDFLIBDIR})
    list(APPEND _LIBS "-L$ENV{NETCDFLIBDIR}")
  endif()
  if (DEFINED ENV{NETCDFLIBNAME})
    list(APPEND _LIBS "$ENV{NETCDFLIBNAME}")
  else()
    # Prefer Fortran lib name
    list(APPEND _LIBS "-lnetcdff" "-lnetcdf")
  endif()
  set(NetCDF_LIBRARIES "${_LIBS}")
  if (NetCDF_INCLUDE_DIRS)
    # Construct imported target from env hints
    add_library(NetCDF::NetCDF_Fortran INTERFACE IMPORTED)
    target_include_directories(NetCDF::NetCDF_Fortran INTERFACE ${NetCDF_INCLUDE_DIRS})
    separate_arguments(NetCDF_LIBRARIES)
    target_link_libraries(NetCDF::NetCDF_Fortran INTERFACE ${NetCDF_LIBRARIES})
    set(_nc_ok TRUE)
    set(NetCDF_Fortran_FOUND TRUE)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF
  REQUIRED_VARS NetCDF_Fortran_FOUND
  FAIL_MESSAGE "NetCDF-Fortran not found. Install libnetcdff-dev (Linux), netcdf-fortran (macOS), or provide hints.")

# Back-compat variables
set(NetCDF_INCLUDE_DIR ${NetCDF_INCLUDE_DIRS})
set(NetCDF_LIBRARY     ${NetCDF_LIBRARIES})
