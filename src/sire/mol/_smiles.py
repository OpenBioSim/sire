__all__ = ["_to_smiles", "_view2d"]

_rdkit_import_error = None

try:
    from rdkit import Chem as _Chem

    _has_rdkit = True
except ImportError:
    _has_rdkit = False
except AttributeError as e:
    _has_rdkit = False
    _rdkit_import_error = e


if _has_rdkit:

    def _to_smiles(obj, include_hydrogens=True, map=None):
        """
        Return this molecule view as a smiles string. Include
        hydrogens in 'include_hydrogens' is True
        """
        from ..convert import sire_to_rdkit
        from ..base import create_map

        map = create_map(map)

        if not include_hydrogens:
            obj = obj[("not element H", map)]

        if not obj.selected_all():
            from ..legacy.Mol import PartialMolecule

            obj = PartialMolecule(obj.molecule(), obj.selection()).extract()

        rdkit_mol = sire_to_rdkit(obj.molecule(), map=map)

        return _Chem.MolToSmiles(rdkit_mol)

    def _view2d(obj, filename=None, map=None):
        """
        Create a 2D representation of this molecule.
        """
        from ..convert import sire_to_rdkit
        from ..base import create_map

        map = create_map(map)

        if not obj.selected_all():
            from ..legacy.Mol import PartialMolecule

            obj = PartialMolecule(obj.molecule(), obj.selection()).extract()

        rdkit_mol = sire_to_rdkit(obj.molecule(), map=map)

        _Chem.AllChem.Compute2DCoords(rdkit_mol)

        from rdkit.Chem import Draw

        if filename is not None:
            Draw.MolToFile(rdkit_mol, filename)

elif _rdkit_import_error is not None:

    def _to_smiles(obj, *args, **kwargs):
        raise ImportError(
            "rdkit cannot be imported. This is because of an error "
            f"when rdkit was loaded ({_rdkit_import_error})."
        )

else:

    def _to_smiles(obj, *args, **kwargs):
        raise ImportError(
            "You need to install rdkit to be able to generate "
            "smiles strings Do this by typing, e.g. "
            "'mamba install -c conda-forge rdkit' and then restarting "
            "Python and running this script/notebook again."
        )
