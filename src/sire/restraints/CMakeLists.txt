########################################
#
# sire.restraints
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        _restraints.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/restraints )
