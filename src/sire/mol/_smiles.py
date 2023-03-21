__all__ = ["_to_smiles", "_view2d", "_selector_to_smiles", "_selector_view2d"]

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

    def _selector_to_smiles(obj, include_hydrogens: bool = False, map=None):
        """
        Return the molecule views in this container as smiles strings. Include
        hydrogens in 'include_hydrogens' is True. This returns a list
        of smiles strings, in the same order as the views in the container
        """
        from ..convert import sire_to_rdkit
        from ..legacy.Convert import rdkit_to_smiles
        from ..base import create_map

        map = create_map(map)

        obj = obj.extract()

        try:
            not_water = obj["not water"].molecules()
        except Exception:
            from . import SelectorMol

            not_water = SelectorMol()

        rdkit_mols = sire_to_rdkit(not_water, map=map)

        if not include_hydrogens:
            from ..legacy.Convert import rdkit_remove_hydrogens

            rdkit_mols = rdkit_remove_hydrogens(rdkit_mols, map)

        smiles = rdkit_to_smiles(rdkit_mols, map)

        if len(not_water) == 1:
            smiles = [smiles]

        try:
            waters = obj["water"].molecules()
        except Exception:
            from . import SelectorMol

            waters = SelectorMol()

        if len(waters) > 0:
            # need to combine them all together
            s = {}

            for mol, smile in zip(not_water, smiles):
                s[mol.number().value()] = smile

            smiles = []

            if include_hydrogens:
                water = "[H]O[H]"
            else:
                water = "O"

            for mol in obj:
                smiles.append(s.get(mol.number().value(), water))

        return smiles

    def _selector_view2d(
        obj,
        filename: str = None,
        height: int = 600,
        width: int = 750,
        include_hydrogens: bool = False,
        num_columns: int = 1,
        force: bool = False,
        map=None,
    ):
        from ..convert import sire_to_rdkit
        from ..legacy.Convert import rdkit_to_smiles
        from ..base import create_map

        if filename is None:
            # we need to be in a jupyter notebook or equivalent
            from IPython.display import SVG

        map = create_map(map)

        obj = obj.extract()

        try:
            not_water = obj["not water"].molecules()
        except Exception:
            from . import SelectorMol

            not_water = SelectorMol()

        # we don't view water molecules
        rdkit_mols = sire_to_rdkit(not_water, map=map)

        # we also don't want any conformers, as these mess up the 2D view
        try:
            for r in rdkit_mols:
                r.RemoveAllConformers()
        except Exception:
            rdkit_mols.RemoveAllConformers()

        if not include_hydrogens:
            from ..legacy.Convert import rdkit_remove_hydrogens

            rdkit_mols = rdkit_remove_hydrogens(rdkit_mols, map)

        smiles = rdkit_to_smiles(rdkit_mols, map)

        # find the unique structures, and count up the rest
        if len(not_water) > 1:
            unique = {}
            unique_mols = []
            ordered_smiles = []

            for smile, rdkit_mol in zip(smiles, rdkit_mols):
                if smile in unique:
                    unique[smile] += 1
                else:
                    unique[smile] = 1
                    unique_mols.append(rdkit_mol)
                    ordered_smiles.append(smile)

            legends = []

            for smile in ordered_smiles:
                legends.append(f"{unique[smile]} x {smile}")
        else:
            unique = {smiles: 1}
            unique_mols = [rdkit_mols]
            legends = [f"1 x {smiles}"]

        # now add back the waters (if any)
        try:
            num_waters = len(obj["water"].molecules())
        except Exception:
            num_waters = 0

        if num_waters > 0:
            from rdkit.Chem import MolFromSmiles

            unique_mols.append(MolFromSmiles("O"))
            unique["O"] = num_waters

            legends.append(f"{num_waters} x O")

        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D

        rdDepictor.SetPreferCoordGen(True)

        num_mols = len(unique_mols)

        if num_columns < 1:
            num_columns = 1

        if num_mols > num_columns:
            num_cols = num_columns
            num_rows = num_mols / num_cols
            if num_mols % num_cols != 0:
                num_rows += 1
        else:
            num_cols = num_mols
            num_rows = 1

        drawer = rdMolDraw2D.MolDraw2DSVG(
            width=width,
            height=height,
            panelWidth=int(width / num_cols),
            panelHeight=int(height / num_rows),
        )

        drawer.DrawMolecules(unique_mols, legends=legends)

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
                        "any other format than SVG. As "
                        "this is not available, we will save the file as a "
                        "SVG. To install `cairosvg` run the command "
                        "'mamba install -c conda-forge cairosvg'"
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

    def _to_smiles(obj, include_hydrogens: bool = False, map=None):
        """
        Return this molecule view as a smiles string. Include
        hydrogens in 'include_hydrogens' is True
        """
        from ..convert import sire_to_rdkit
        from ..legacy.Convert import rdkit_to_smiles
        from ..base import create_map

        # check if this is a water molecule
        if obj.selected_all():
            try:
                obj = obj["water"]
                if obj.selected_all():
                    if include_hydrogens:
                        return "[H]O[H]"
                    else:
                        return "O"
            except Exception:
                pass

        map = create_map(map)

        rdkit_mol = sire_to_rdkit(obj.extract(), map)

        if not include_hydrogens:
            from ..legacy.Convert import rdkit_remove_hydrogens

            rdkit_mol = rdkit_remove_hydrogens(rdkit_mol, map)

        return rdkit_to_smiles(rdkit_mol, map)

    def _view2d(
        obj,
        filename: str = None,
        height: int = 300,
        width: int = 900,
        include_hydrogens: bool = False,
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

        rdkit_mol = None

        # check for water
        if obj.selected_all():
            try:
                obj = obj["water"]
                if obj.selected_all():
                    from ..legacy.Convert import smiles_to_rdkit

                    rdkit_mol = smiles_to_rdkit("O", "O", map)
            except Exception:
                pass

        if rdkit_mol is None:
            rdkit_mol = sire_to_rdkit(obj.extract(), map)
            # remove conformers, as these mess up the 2D view
            rdkit_mol.RemoveAllConformers()

        if not include_hydrogens:
            from ..legacy.Convert import rdkit_remove_hydrogens

            rdkit_mol = rdkit_remove_hydrogens(rdkit_mol, map)

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
                        "any other format than SVG. As "
                        "this is not available, we will save the file as a "
                        "SVG. To install `cairosvg` run the command "
                        "'mamba install -c conda-forge cairosvg'"
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

    def _no_rdkit():
        raise ImportError(
            "rdkit cannot be imported. This is because of an error "
            f"when rdkit was loaded ({_rdkit_import_error})."
        )

    def _view2d(obj, *args, **kwargs):
        _no_rdkit()

    def _to_smiles(obj, *args, **kwargs):
        _no_rdkit()

    def _selector_view2d(obj, *args, **kwargs):
        _no_rdkit()

    def _selector_to_smiles(obj, *args, **kwargs):
        _no_rdkit()

else:

    def _view2d(obj, *args, **kwargs):
        raise ImportError(
            "You need to install rdkit to be able to generate "
            "2D views of molecules. Do this by typing, e.g. "
            "'mamba install -c conda-forge rdkit' and then restarting "
            "Python and running this script/notebook again."
        )

    def _to_smiles(obj, *args, **kwargs):
        raise ImportError(
            "You need to install rdkit to be able to generate "
            "smiles strings Do this by typing, e.g. "
            "'mamba install -c conda-forge rdkit' and then restarting "
            "Python and running this script/notebook again."
        )

    def _selector_to_smiles(obj, *args, **kwargs):
        _to_smiles(obj, *args, **kwargs)

    def _selector_view2d(obj, *args, **kwargs):
        _view2d(obj, *args, **kwargs)
