########################################
#
# sire.io
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/io )

add_subdirectory (parser)
