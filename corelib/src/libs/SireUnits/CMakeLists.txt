########################################
#
# CMake file for library:SireUnits
#
########################################

# Export symbols when building the library
add_definitions( "-DSIREUNITS_BUILD" )

# Third Party dependencies of this module

# Other Sire libraries
include_directories( ${CMAKE_SOURCE_DIR}/src/libs )

set ( SIREUNITS_HEADERS
      angle.h
      ast.h
      convert.h
      dimensions.h
      generalunit.h
      temperature.h
      units.h
    )

# Define the sources in SireUnits
set ( SIREUNITS_SOURCES

      register_sireunits.cpp

      ast.cpp
      angle.cpp
      convert.cpp
      dimensions.cpp
      generalunit.cpp
      parser.cpp
      temperature.cpp

      ${SIREUNITS_HEADERS}
    )

add_library (SireUnits ${SIREUNITS_SOURCES})

set_target_properties (SireUnits PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireUnits
                       SireStream
                       )

# installation
install( TARGETS SireUnits EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREUNITS_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireUnits )
