__all__ = [
    "sire_to_rdkit",
    "rdkit_to_sire",
    "rdkit_to_smiles",
    "rdkit_remove_hydrogens",
    "smiles_to_rdkit",
    "sire_to_openmm",
    "openmm_to_sire",
    "supported_formats",
]

try:
    from ._SireRDKit import (
        sire_to_rdkit,
        rdkit_to_sire,
        rdkit_to_smiles,
        rdkit_remove_hydrogens,
        smiles_to_rdkit,
    )

    _has_rdkit = True

except Exception:
    # RDKit support is not available
    def _no_rdkit():
        raise ModuleNotFoundError(
            "Unable to convert to/from RDKit as it is not installed. "
            "Please install using `mamba install -c conda-forge rdkit` "
            "and then re-run this script."
        )

    _has_rdkit = False

    def sire_to_rdkit(*args, **kwargs):
        _no_rdkit()

    def rdkit_to_sire(*args, **kwargs):
        _no_rdkit()

    def rdkit_to_smiles(*args, **kwargs):
        _no_rdkit()

    def rdkit_remove_hydrogens(*args, **kwargs):
        _no_rdkit()

    def smiles_to_rdkit(*args, **kwargs):
        _no_rdkit()


try:
    from ._SireOpenMM import (
        _sire_to_openmm_system,
        _openmm_system_to_sire,
    )

    _has_openmm = True

    def openmm_to_sire(mols, map):
        import openmm

        if type(mols) is not openmm.System:
            raise TypeError(
                "You can only convert an openmm.System to sire, not "
                f"a {type(mols)}."
            )

        # Need to be sure that 'mols' is an openmm.System or else
        # we will crash!
        return _openmm_system_to_sire(mols, map)

    def sire_to_openmm(mols, map):
        import openmm

        system = openmm.System()

        # system must be an openmm.System() or else we will crash!
        _sire_to_openmm_system(system, mols, map)

        return system

except Exception:
    # OpenMM support is not available
    def _no_openmm():
        raise ModuleNotFoundError(
            "Unable to convert to/from OpenMM as this code hasn't been "
            "written yet. We hope to support this soon!"
        )

    _has_openmm = False

    def sire_to_openmm(*args, **kwargs):
        _no_openmm()

    def openmm_to_sire(*args, **kwargs):
        _no_openmm()


def supported_formats():
    """Return all of the formats supported by this installation"""
    f = ["sire"]

    if _has_openmm:
        f.append("openmm")

    if _has_rdkit:
        f.append("rdkit")

    import sys

    # BioSimSpace needs to have already been loaded
    # otherwise it can't be imported
    if "BioSimSpace" in sys.modules:
        f.append("biosimspace")

    f.sort()

    return f
