########################################
#
# CMake file for Converters
#
########################################

# Set the compiler flags so that we compile optimised, rather than
# small code
set( CMAKE_CXX_FLAGS "${OPTIMISED_CXX_FLAGS}" )

add_subdirectory (SireRDKit)
add_subdirectory (SireOpenMM)
#add_subdirectory (SireGemmi)

# installation
set( INSTALLDIR ${SIRE_PYTHON}/sire/legacy/Convert )

install( FILES __init__.py
         DESTINATION ${INSTALLDIR}
       )

