########################################
#
# CMake file for Converters
#
########################################

add_subdirectory (SireRDKit)

# installation
set( INSTALLDIR ${SIRE_PYTHON}/sire/legacy/Convert )

install( FILES __init__.py
         DESTINATION ${INSTALLDIR}
       )

