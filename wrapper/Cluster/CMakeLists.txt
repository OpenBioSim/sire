########################################
#
# CMake file for Python wrapper:_Cluster
#
########################################

# Third Party dependencies of this module
include_directories( ${BOOST_INCLUDE_DIR}
                     ${PYTHON_INCLUDE_DIR} )

if ( SIRE_USE_MPI )
  include_directories( BEFORE SYSTEM ${SIRE_MPI_INCLUDE_PATH} )
  add_definitions( ${SIRE_MPI_DEFINITIONS} )
endif()

# Sire include paths
include_directories( BEFORE ${SIRE_INCLUDE_DIR}/SireCluster
                            ${SIRE_INCLUDE_DIR} )

# Other python wrapping directories
include_directories(${CMAKE_SOURCE_DIR})

# Get the list of autogenerated files
include(CMakeAutogenFile.txt)

# Define the sources in Cluster
set ( WRAPCLUSTER_SOURCES

      _Cluster.main.cpp

      ${PYPP_SOURCES}

    )

# Create the library that holds all of the class wrappers
add_library (Cluster ${WRAPCLUSTER_SOURCES})

set_target_properties (Cluster
                       PROPERTIES
                         VERSION ${SIRE_VERSION}
                         SOVERSION ${SIRE_VERSION_MAJOR}
                         OUTPUT_NAME "Cluster"
                         PREFIX "_"
                         SUFFIX ${PYTHON_MODULE_EXTENSION}
		       )

target_link_libraries (Cluster
                       SirePython
                       SIRE_SireCluster
                       ${Boost_LIBRARIES}
                      )

include( LimitSirePythonExportSymbols )
export_this_symbol_only( "init_Cluster" "_init_Cluster" )

# installation
set( INSTALLDIR ${SIRE_PYTHON}/sire/legacy/Cluster )

install( TARGETS Cluster
         LIBRARY DESTINATION ${INSTALLDIR}
         ARCHIVE DESTINATION ${INSTALLDIR}
         RUNTIME DESTINATION ${INSTALLDIR}
        )

install( FILES __init__.py
         DESTINATION ${INSTALLDIR}
       )

