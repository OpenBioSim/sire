########################################
#
# CMake file for library:SireIO
#
########################################

# Export symbols when building the library
add_definitions( "-DSIREIO_BUILD" )

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# This library uses Intel Threaded Building blocks
include_directories(${TBB_INCLUDE_DIR})

# If we have NetCDF, then include the header files
if ( ${SIRE_FOUND_NETCDF} )
  message( STATUS "Adding in NetCDF support to SireIO using ${NetCDF_ROOT_DIR}" )
  message( STATUS "${NetCDF_INCLUDE_DIR} | ${NetCDF_LIBRARIES}" )
  include_directories(${NetCDF_INCLUDE_DIR})
  set ( SIREIO_NETCDF_LIBRARIES ${NetCDF_LIBRARIES} )
  add_definitions( "-DSIRE_USE_NETCDF" )
endif()

# Define the headers in SireIO
set ( SIREIO_HEADERS
      amber.h
      amberformat.h
      amberprm.h
      amberrst.h
      amberrst7.h
      ambertraj.h
      biosimspace.h
      charmmpsf.h
      cube.h
      dcd.h
      errors.h
      filetrajectory.h
      filetrajectoryparser.h
      flexibilitylibrary.h
      fortranfile.h
      getcoordsarray.h
      gro87.h
      grotop.h
      iobase.h
      mol2.h
      moleculeparser.h
      netcdffile.h
      pdb.h
      pdb2.h
      pdbx.h
      perturbationslibrary.h
      protoms.h
      sdf.h
      supplementary.h
      textfile.h
      tinker.h
      trajectorymonitor.h
      trr.h
      xdrfile.h
      xtc.h
      zmatrixmaker.h
    )

# Define the sources in SireIO
set ( SIREIO_SOURCES

      register_sireio.cpp

      third_party/xdrfile_trr.c
      third_party/xdrfile_xtc.c
      third_party/xdrfile.c
      third_party/xtc_seek.c

      amber.cpp
      amberformat.cpp
      amberprm.cpp
      amberrst.cpp
      amberrst7.cpp
      ambertraj.cpp
      biosimspace.cpp
      charmmpsf.cpp
      cube.cpp
      dcd.cpp
      errors.cpp
      filetrajectory.cpp
      filetrajectoryparser.cpp
      flexibilitylibrary.cpp
      fortranfile.cpp
      getcoordsarray.cpp
      gro87.cpp
      grotop.cpp
      iobase.cpp
      mol2.cpp
      moleculeparser.cpp
      netcdffile.cpp
      pdb.cpp
      pdb2.cpp
      pdbx.cpp
      perturbationslibrary.cpp
      protoms.cpp
      sdf.cpp
      supplementary.cpp
      textfile.cpp
      tinker.cpp
      trajectorymonitor.cpp
      trr.cpp
      xdrfile.cpp
      xtc.cpp
      zmatrixmaker.cpp

      ${SIREIO_HEADERS}
    )

add_library (SireIO ${SIREIO_SOURCES})

set_target_properties (SireIO PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireIO
                       SireMM
                       SireMove
                       SireMol
                       SireUnits
                       SireStream
                       ${SIREIO_NETCDF_LIBRARIES}
                       )

# installation
install( TARGETS SireIO EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREIO_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireIO )
