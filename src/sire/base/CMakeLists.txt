########################################
#
# sire.base
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        _progressbar.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/base )
