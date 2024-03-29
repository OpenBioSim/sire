########################################
#
# CMake file for library:SireMaths
#
########################################

# Export symbols when building the library
add_definitions( "-DSIREMATHS_BUILD" )

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# Define the headers in SireMaths
set ( SIREMATHS_HEADERS
      accumulator.h
      align.h
      axisset.h
      boys.h
      complex.h
      constants.h
      distvector.h
      energytrajectory.h
      errors.h
      freeenergyaverage.h
      gamma.h
      histogram.h
      line.h
      linearap.h
      maths.h
      matrix.h
      multidouble.h
      multifixed.h
      multifloat.h
      multiint.h
      multiquaternion.h
      multiuint.h
      multivector.h
      nmatrix.h
      nvector.h
      n4matrix.h
      plane.h
      quaternion.h
      rangenerator.h
      rational.h
      rotate.h
      sincos.h
      sphere.h
      sphereproperty.h
      torsion.h
      triangle.h
      trigmatrix.h
      vector.h
      vector3d.hpp
      vectorproperty.h

      third_party/eig3/eig3.cpp
    )

# Define the sources in SireMaths
set ( SIREMATHS_SOURCES

      register_siremaths.cpp

      accumulator.cpp
      align.cpp
      axisset.cpp
      boys.cpp
      complex.cpp
      distvector.cpp
      energytrajectory.cpp
      errors.cpp
      freeenergyaverage.cpp
      gamma.cpp
      histogram.cpp
      line.cpp
      linearap.cpp
      matrix.cpp
      multidouble.cpp
      multifixed.cpp
      multifloat.cpp
      multiint.cpp
      multiquaternion.cpp
      multiuint.cpp
      multivector.cpp
      nmatrix.cpp
      nvector.cpp
      n4matrix.cpp
      plane.cpp
      quaternion.cpp
      rangenerator.cpp
      rational.cpp
      sincos.cpp
      sphere.cpp
      sphereproperty.cpp
      torsion.cpp
      triangle.cpp
      trigmatrix.cpp
      vector.cpp
      vectorproperty.cpp

      ${SIREMATHS_HEADERS}
    )

add_library (SireMaths ${SIREMATHS_SOURCES})

set_target_properties (SireMaths PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries, and GSL
target_link_libraries (SireMaths
                       SireID
                       SireBase
                       SireUnits
                       SireStream
                       ${SIREMATHS_EXTRA_LIBRARIES}
                      )

save_sire_variable( "SIREMATHS_EXTRA_LIBRARIES" "${SIREMATHS_EXTRA_LIBRARIES}" )

# installation
install( TARGETS SireMaths EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREMATHS_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireMaths )
