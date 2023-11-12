##########################################
#
# CMake file for Python wrapper:_SireGemmi
#
##########################################

# Find gemmi and pybind11 libraries
find_package(gemmi)
find_package(pybind11)

if (${gemmi_FOUND} AND ${pybind11_FOUND})
  # Only compile and install if we have gemmi
  message( STATUS "Compiling SireGemmi converter" )

  get_target_property(GEMMI_LOCATION gemmi::gemmi_cpp LOCATION)

  message( STATUS "GEMMI LIBRARY ${GEMMI_LOCATION}" )
  message( STATUS "PYBIND11 INCLUDE ${pybind11_INCLUDE_DIR}")
  message( STATUS "PYBIND11 LIBRARY ${pybind11_LIBRARIES}" )

  # Third Party dependencies of this module
  include_directories( ${PYTHON_INCLUDE_DIR}
                       ${pybind11_INCLUDE_DIR}
                       ${gemmi_INCLUDE_DIR} )

  # Sire include paths
  include_directories( BEFORE ${SIRE_INCLUDE_DIR} )

  # Other python wrapping directories
  include_directories(${CMAKE_SOURCE_DIR})

  # Define the sources in SireRDKit
  set ( SIREGEMMI_SOURCES

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
                          gemmi::gemmi_cpp
                          ${pybind11_LIBRARIES}
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
