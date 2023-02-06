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

    def _view2d(
        obj,
        filename: str = None,
        height: int = 300,
        width: int = 900,
        map=None,
    ):
        """
        Create a 2D representation of this molecule. If 'filename'
        is set then this will be written to that file. Otherwise
        this will be returned for visualisation in a jupyter notebook.
        """
        from ..convert import sire_to_rdkit
        from ..base import create_map

        if filename is None:
            # we need to be in a jupyter notebook or equivalent
            from IPython.display import SVG

        map = create_map(map)

        if not obj.selected_all():
            from ..legacy.Mol import PartialMolecule

            obj = PartialMolecule(obj.molecule(), obj.selection()).extract()

        rdkit_mol = sire_to_rdkit(obj.molecule(), map=map)

        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D

        rdDepictor.SetPreferCoordGen(True)

        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(rdkit_mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        if filename is not None:
            import os

            (basename, format) = os.path.splitext(os.path.abspath(filename))

            while format.startswith("."):
                format = format[1:]

            if len(format) == 0:
                format = "svg"

            format = format.lower()

            if format != "svg":
                try:
                    import cairosvg
                except Exception:
                    from ..utils import Console

                    Console.warning(
                        "We need the module `cairosvg` to save image files in "
                        "any other format than SVG. . As "
                        "this is not available, we will save the file as a SVG"
                    )

                    format = "svg"

            filename = f"{basename}.{format}"

            if format == "svg":
                with open(filename, "w") as FILE:
                    FILE.write(svg)

            elif format == "pdf":
                cairosvg.svg2pdf(bytestring=svg.encode(), write_to=filename)

            elif format == "png":
                cairosvg.svg2png(
                    bytestring=svg.encode(),
                    write_to=filename,
                    output_width=width,
                    output_height=height,
                )

            else:
                raise NotImplementedError(
                    f"Cannot save image as format {format}. This is not yet "
                    "supported in this code. Please choose 'png', 'pdf' or "
                    "'svg'."
                )

            return filename

        else:
            return SVG(svg)

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
