########################################
#
# CMake file for Python wrapper:_Mol
#
########################################

# Third Party dependencies of this module
include_directories( ${BOOST_INCLUDE_DIR}
                     ${PYTHON_INCLUDE_DIR}
                     ${SIRE_GSL_INCLUDE_PATH} )

# Sire include paths
include_directories( BEFORE ${SIRE_INCLUDE_DIR}/SireMol
                            ${SIRE_INCLUDE_DIR}/SireMM
                            ${SIRE_INCLUDE_DIR} )

# Other python wrapping directories
include_directories(${CMAKE_SOURCE_DIR})

# Get the list of autogenerated files
include(CMakeAutogenFile.txt)

# Define the sources in Mol
set ( WRAPMOL_SOURCES

      _Mol.main.cpp

      ${PYPP_SOURCES}

    )

# Create the library that holds all of the class wrappers
add_library (Mol ${WRAPMOL_SOURCES})

set_target_properties (Mol
                       PROPERTIES
                         VERSION ${SIRE_VERSION}
                         SOVERSION ${SIRE_VERSION_MAJOR}
                         OUTPUT_NAME "Mol"
                         PREFIX "_"
                         SUFFIX ${PYTHON_MODULE_EXTENSION}
		       )

target_link_libraries (Mol
                       SirePython
                       SIRE_SireID
                       SIRE_SireMol
                       SIRE_SireMM
                       ${Boost_LIBRARIES}
                      )

include( LimitSirePythonExportSymbols )
export_this_symbol_only( "init_Mol" "_init_Mol" )

# installation
set( INSTALLDIR ${SIRE_PYTHON}/sire/legacy/Mol )

install( TARGETS Mol
         LIBRARY DESTINATION ${INSTALLDIR}
         ARCHIVE DESTINATION ${INSTALLDIR}
         RUNTIME DESTINATION ${INSTALLDIR}
        )

install( FILES __init__.py
         DESTINATION ${INSTALLDIR}
       )

