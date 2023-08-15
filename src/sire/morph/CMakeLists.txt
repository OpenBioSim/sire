########################################
#
# sire.morph
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        _ghost_atoms.py
        _perturbation.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/morph )
