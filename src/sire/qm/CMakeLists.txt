########################################
#
# sire.qm
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        _emle.py
        _utils.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/qm )
