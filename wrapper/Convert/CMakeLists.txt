########################################
#
# CMake file for Converters
#
########################################

add_subdirectory (SireRDKit)
add_subdirectory (SireOpenMM)

# installation
set( INSTALLDIR ${SIRE_PYTHON}/sire/legacy/Convert )

install( FILES __init__.py
         DESTINATION ${INSTALLDIR}
       )

