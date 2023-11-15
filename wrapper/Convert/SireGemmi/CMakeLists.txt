##########################################
#
# CMake file for Python wrapper:_SireGemmi
#
##########################################

# Find gemmi and pybind11 libraries
find_package(gemmi)
find_package(pybind11)

# gemmi needs zlib, at least on Linux. Make sure we find and
# link to this
cmake_policy(SET CMP0074 NEW)
set( ZLIB_ROOT "${ANACONDA_BASE}" NO_DEFAULT_PATH)
find_package(ZLIB)

if (${ZLIB_FOUND})
  # test if ${ZLIB_LIBRARIES} starts with ${ANACONDA_BASE}
  if (NOT ${ZLIB_LIBRARIES} MATCHES "^${ANACONDA_BASE}")
    message( STATUS "ZLIB found, but not in ${ANACONDA_BASE}. This may cause problems." )
    set( ZLIB_LIBRARIES "${ANACONDA_BASE}/lib/libz.a" )
    if (EXISTS "${ZLIB_LIBRARIES}")
      message("Manually using ${ZLIB_LIBRARIES} for ZLIB_LIBRARIES")
    else()
      message( WARNING "ZLIB not found in ${ANACONDA_BASE}/lib." )
      set(ZLIB_FOUND FALSE)
      unset(ZLIB_LIBRARIES)
    endif()
  endif()
endif()

if (${gemmi_FOUND} AND ${pybind11_FOUND})
  # Only compile and install if we have gemmi
  message( STATUS "Compiling SireGemmi converter" )

  get_target_property(GEMMI_LIBRARY gemmi::gemmi_cpp LOCATION)

  message( STATUS "GEMMI LIBRARY ${GEMMI_LIBRARY}" )
  message( STATUS "PYBIND11 INCLUDE ${pybind11_INCLUDE_DIR}")
  message( STATUS "PYBIND11 LIBRARY ${pybind11_LIBRARIES}" )
  message( STATUS "ZLIB LIBRARIES ${ZLIB_LIBRARIES}" )

  # Third Party dependencies of this module
  include_directories( ${PYTHON_INCLUDE_DIR}
                       ${pybind11_INCLUDE_DIR}
                       ${gemmi_INCLUDE_DIR} )

  # Sire include paths
  include_directories( BEFORE ${SIRE_INCLUDE_DIR} )

  # Other python wrapping directories
  include_directories(${CMAKE_SOURCE_DIR})

  # Do we need to manually add in the functions from GEMMI_EXTRA_SOURCE?
  set(OLD_CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
  set(CMAKE_REQUIRED_LIBRARIES ${GEMMI_LIBRARY})
  check_cxx_source_runs(
    "
    #include <gemmi/cif.hpp>
    #include <gemmi/mmcif.hpp>
    int main() {
      gemmi::CifFile cf;
      gemmi::MmcifFile mf;
      return 0;
    }
    "
    GEMMI_EXTRA_SOURCE_WORKS
  )
  set(CMAKE_REQUIRED_LIBRARIES ${OLD_CMAKE_REQUIRED_LIBRARIES})

  if (GEMMI_EXTRA_SOURCE_WORKS)
    message( STATUS "Gemmi library links ok and contains the necessary functions.")
    unset( GEMMI_EXTRA_SOURCE )
  else()
    message( STATUS "GEMMI_EXTRA_SOURCE_WORKS is false, so we need to add in the extra source files" )
    set ( GEMMI_EXTRA_SOURCE
          # need this as it isn't linked in property from the gemmi static lib?
          third_party/mmcif.cpp
          third_party/to_mmcif.cpp
          third_party/sprintf.cpp
          third_party/polyheur.cpp
          third_party/resinfo.cpp
        )
  endif()

  # Define the sources in SireRDKit
  set ( SIREGEMMI_SOURCES

        ${GEMMI_EXTRA_SOURCE}

        _SireGemmi.main.cpp

        sire_gemmi.cpp

      )

  # Create the library that holds all of the class wrappers
  add_library (SireGemmi ${SIREGEMMI_SOURCES})

  set_target_properties(SireGemmi
                        PROPERTIES
                          VERSION ${SIRE_VERSION}
                          SOVERSION ${SIRE_VERSION_MAJOR}
                          OUTPUT_NAME "SireGemmi"
                          PREFIX "_"
                          SUFFIX ${PYTHON_MODULE_EXTENSION}
                        )

  target_link_libraries(SireGemmi
                          SirePython
                          SIRE_SireIO
                          SIRE_SireMM
                          SIRE_SireMol
                          SIRE_SireBase
                          SIRE_SireStream
                          SIRE_SireError
                          ${GEMMI_LIBRARY}
                          ${pybind11_LIBRARIES}
                          ${ZLIB_LIBRARIES}
                        )


  include( LimitSirePythonExportSymbols )
  export_this_symbol_only( "init_SireRDKit" "_init_SireGemmi" )

  # installation
  set( INSTALLDIR ${SIRE_PYTHON}/sire/legacy/Convert )

  install( TARGETS SireGemmi
          LIBRARY DESTINATION ${INSTALLDIR}
          ARCHIVE DESTINATION ${INSTALLDIR}
          RUNTIME DESTINATION ${INSTALLDIR}
          )
else()
  message( WARNING "gemmi not found, so we cannot compile the converter." )
endif()
