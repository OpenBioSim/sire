##########################################
#
# CMake file for Python wrapper:_SireRDKit
#
##########################################

# Find RDKit libraries
find_package(RDKit REQUIRED)

if (${RDKIT_FOUND})
  # Only compile and install if we have RDKit
  message( STATUS "Compiling SireRDKit converter" )

  message( STATUS "RDKIT_INCLUDE_DIR: ${RDKIT_INCLUDE_DIR}" )
  message( STATUS "RDKIT_LIBRARIES:   ${RDKIT_LIBRARIES}" )

  # Third Party dependencies of this module
  include_directories( ${BOOST_INCLUDE_DIR}
                       ${PYTHON_INCLUDE_DIR}
                       ${RDKIT_INCLUDE_DIR} )

  # Sire include paths
  include_directories( BEFORE ${SIRE_INCLUDE_DIR} )

  # Other python wrapping directories
  include_directories(${CMAKE_SOURCE_DIR})

  # Define the sources in SireRDKit
  set ( SIRERDKIT_SOURCES

        _SireRDKit.main.cpp

        sire_rdkit.cpp

      )

  # Create the library that holds all of the class wrappers
  add_library (SireRDKit ${SIRERDKIT_SOURCES})

  set_target_properties(SireRDKit
                         PROPERTIES
                           VERSION ${SIRE_VERSION}
                           SOVERSION ${SIRE_VERSION_MAJOR}
                           OUTPUT_NAME "SireRDKit"
                           CXX_STANDARD 17
                           PREFIX "_"
                           SUFFIX ${PYTHON_MODULE_EXTENSION}
                        )

  target_link_libraries(SireRDKit
                          SirePython
                          SIRE_SireMol
                          SIRE_SireBase
                          SIRE_SireStream
                          SIRE_SireError
                          ${Boost_LIBRARIES}
                          ${RDKIT_LIBRARIES}
                        )


  include( LimitSirePythonExportSymbols )
  export_this_symbol_only( "init_SireRDKit" "_init_SireRDKit" )

  # installation
  set( INSTALLDIR ${SIRE_PYTHON}/sire/legacy/Convert )

  install( TARGETS SireRDKit
          LIBRARY DESTINATION ${INSTALLDIR}
          ARCHIVE DESTINATION ${INSTALLDIR}
          RUNTIME DESTINATION ${INSTALLDIR}
          )
else()
  message( WARNING "RDKit not found, so we cannot compile the converter." )
endif()
