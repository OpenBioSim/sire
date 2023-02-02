__all__ = ["sire_to_rdkit", "rdkit_to_sire"]


try:
    from ._SireRDKit import sire_to_rdkit, rdkit_to_sire
except Exception:
    # RDKit support is not available
    def _no_rdkit():
        raise ModuleNotFoundError(
            "Unable to convert to/from RDKit as it is not installed. "
            "Please install using `mamba install -c conda-forge rdkit` "
            "and then re-run this script."
        )

    def sire_to_rdkit(*args, **kwargs):
        _no_rdkit()

    def rdkit_to_sire(*args, **kwargs):
        _no_rdkit()
