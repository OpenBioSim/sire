##########################################
#
# CMake file for Python wrapper:_SireRDKit
#
##########################################

if (SIRE_HAS_CPP_17)
  # Find RDKit libraries
  find_package(RDKit)

  if (${RDKIT_FOUND})
    # Only compile and install if we have RDKit
    message( STATUS "Compiling SireRDKit converter" )

    message( STATUS "RDKIT_INCLUDE_DIR: ${RDKIT_INCLUDE_DIR}" )
    message( STATUS "RDKIT_LIBRARIES:   ${RDKIT_LIBRARIES}" )

    # Third Party dependencies of this module
    include_directories( ${BOOST_INCLUDE_DIR}
                         ${PYTHON_INCLUDE_DIR}
                         ${RDKIT_INCLUDE_DIR} )

    if (WIN32)
      # RDKit needs this definition set on Windows
      add_definitions(-DWIN32)
      add_definitions(-DRDKIT_DYN_LINK)
    endif()

    if (APPLE)
      # RDKit uses std::bad_any_cast, which isn't available by default
      # with MacOS 10.9. But it is supplied with conda, so need to mark
      # that, e.g. see
      # https://conda-forge.org/docs/maintainer/knowledge_base.html#newer-c-features-with-old-sdk
      add_definitions(-D_LIBCPP_DISABLE_AVAILABILITY)
    endif()

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
                            PREFIX "_"
                            SUFFIX ${PYTHON_MODULE_EXTENSION}
                          )

    target_link_libraries(SireRDKit
                            SirePython
                            SIRE_SireMol
                            SIRE_SireMM
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
else()
  message( WARNING "RDKit support disabled as the compiler doesn't support C++ 2017")
endif()
