########################################
#
# sire.morph
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        _alchemy.py
        _ghost_atoms.py
        _hmr.py
        _pertfile.py
        _perturbation.py
        _repex.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/morph )
