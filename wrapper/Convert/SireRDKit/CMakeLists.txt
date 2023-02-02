##########################################
#
# CMake file for Python wrapper:_SireRDKit
#
##########################################

# Third Party dependencies of this module
include_directories( ${BOOST_INCLUDE_DIR}
                     ${PYTHON_INCLUDE_DIR}
                     ${RDKIT_INCLUDE_PATH} )

# Sire include paths
include_directories( BEFORE ${SIRE_INCLUDE_DIR} )

# Other python wrapping directories
include_directories(${CMAKE_SOURCE_DIR})

# Define the sources in SireRDKit
set ( SIRERDKIT_SOURCES

      _SireRDKit.main.cpp

    )

# Create the library that holds all of the class wrappers
add_library (SireRDKit ${SIRERDKIT_SOURCES})

set_target_properties (SireRDKit
                       PROPERTIES
                         VERSION ${SIRE_VERSION}
                         SOVERSION ${SIRE_VERSION_MAJOR}
                         OUTPUT_NAME "SireRDKit"
                         PREFIX "_"
                         SUFFIX ${PYTHON_MODULE_EXTENSION}
		       )

target_link_libraries (SireRDKit
                       SirePython
                       SIRE_SireMol
                       ${Boost_LIBRARIES}
                       ${Rdkit_LIBRARIES}
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
