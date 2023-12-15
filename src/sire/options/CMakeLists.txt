########################################
#
# sire.options
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        _option.py
        _dynamics_options.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/options )
