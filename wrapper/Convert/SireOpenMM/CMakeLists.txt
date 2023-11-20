###########################################
#
# CMake file for Python wrapper:_SireOpenMM
#
###########################################

if (${SIRE_USE_OPENMM})
  # look for OpenMM in the system
  set(OpenMM_ROOT_DIR "${ANACONDA_BASE}")
  set(OPENMM_ROOT_DIR "${ANACONDA_BASE}")
  FIND_PACKAGE( OpenMM )

  if ( ${OpenMM_FOUND} )
    message ( STATUS "Adding in OpenMM support to SireConvert using OpenMM ${OpenMM_ROOT_DIR}" )
    message ( STATUS "${OpenMM_LIBRARIES} | ${OpenMM_INCLUDE_DIR}" )

    # Third Party dependencies of this module
    include_directories( ${OpenMM_INCLUDE_DIR} )
    add_definitions("-DSIRE_USE_OPENMM")
    add_definitions(-DQT_NO_SIGNALS_SLOTS_KEYWORDS)

    # Sire include paths
    include_directories( BEFORE ${SIRE_INCLUDE_DIR} )

    # Other python wrapping directories
    include_directories(${CMAKE_SOURCE_DIR})

    include(CheckCXXSourceRuns)

    # Whether or not we have CustomCPPForceImpl (OpenMM 8.1+)
    CHECK_CXX_SOURCE_RUNS("
          #include \"openmm/internal/ContextImpl.h\"
          #include \"openmm/internal/CustomCPPForceImpl.h\"

          int main()
          {
              OpenMM::ContextImpl* context = NULL;
              OpenMM::CustomCPPForceImpl* force = NULL;
              return 0;
          }"
        CAN_USE_CUSTOMCPPFORCE)

    if (CAN_USE_CUSTOMCPPFORCE)
      message(STATUS "OpenMM version supports CustomCPPForce")
      add_definitions("-DSIRE_USE_CUSTOMCPPFORCE")
    else()
      message(STATUS "OpenMM version does not support CustomCPPForce")
    endif()

    # Define the sources in SireOpenMM
    set ( SIREOPENMM_SOURCES

          _SireOpenMM.main.cpp

          customforce.cpp
          emle.cpp
          lambdalever.cpp
          openmmminimise.cpp
          openmmmolecule.cpp
          qmmm.cpp
          sire_openmm.cpp
          sire_to_openmm_system.cpp

          lbgfs/lbfgs.cpp

      )

    # Create the library that holds all of the class wrappers
    add_library (SireOpenMM ${SIREOPENMM_SOURCES})

    set_target_properties(SireOpenMM
                          PROPERTIES
                          VERSION ${SIRE_VERSION}
                          SOVERSION ${SIRE_VERSION_MAJOR}
                          OUTPUT_NAME "SireOpenMM"
                          PREFIX "_"
                          SUFFIX ${PYTHON_MODULE_EXTENSION}
                          )

    target_link_libraries(SireOpenMM
                          SirePython
                          SIRE_SireMol
                          SIRE_SireMM
                          SIRE_SireSystem
                          SIRE_SireCAS
                          SIRE_SireBase
                          SIRE_SireStream
                          SIRE_SireError
                          ${OpenMM_LIBRARIES}
                          )

    include( LimitSirePythonExportSymbols )
    export_this_symbol_only( "init_SireOpenMM" "_init_SireOpenMM" )

    # installation
    set( INSTALLDIR ${SIRE_PYTHON}/sire/legacy/Convert )

    install( TARGETS SireOpenMM
          LIBRARY DESTINATION ${INSTALLDIR}
          ARCHIVE DESTINATION ${INSTALLDIR}
          RUNTIME DESTINATION ${INSTALLDIR}
          )

    install( FILES _sommcontext.py
         DESTINATION ${INSTALLDIR}
       )
  else()
    message ( STATUS "OpenMM not found, so OpenMM support is disabled." )
  endif()
else()
  message( STATUS "OpenMM support disabled by user" )
endif()
