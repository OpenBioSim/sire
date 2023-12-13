##########################################
#
# CMake file for Python wrapper:_SireGemmi
#
##########################################

# Find gemmi and pybind11 libraries
find_package(gemmi)
find_package(pybind11)

if (${gemmi_FOUND})
  # We can only link against the dynamic library - static doesn't work
  get_target_property(GEMMI_LIBRARY gemmi::gemmi_cpp LOCATION)

  # Test if we can link to Gemmi on Windows
  if (WIN32)
    # test if GEMMI_LIBRARY ends with .lib
    if (${GEMMI_LIBRARY} MATCHES "\\.lib$")
      message( STATUS "GEMMI_LIBRARY is a static library." )
      message( STATUS "We can only link against the dynamic library." )
      message( STATUS "Disabling gemmi PDBx/mmCIF support.")
      set(gemmi_FOUND FALSE)
    endif()
  else()
    # test if GEMMI_LIBRARY ends with .a
    if (${GEMMI_LIBRARY} MATCHES "\\.a$")
      message( STATUS "GEMMI_LIBRARY is a static library." )
      message( STATUS "We can only link against the dynamic library." )
      message( STATUS "Disabling gemmi PDBx/mmCIF support.")
      set(gemmi_FOUND FALSE)
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
                          ${GEMMI_LIBRARY}
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
