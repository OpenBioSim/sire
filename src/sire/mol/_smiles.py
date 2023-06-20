__all__ = [
    "_to_smiles",
    "_to_smarts",
    "_view2d",
    "_selector_to_smiles",
    "_selector_to_smarts",
    "_selector_view2d",
]

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

    def _selector_to_smarts(
        obj, as_search: bool = False, include_hydrogens: bool = False, map=None
    ):
        """
        Return the molecule views in this container as smarts strings. Include
        hydrogens if 'include_hydrogens' is True. This returns a list of
        smarts strings, in the same order as the views in the container.

        If 'as_search' is True then the smarts string is
        returned as a sire search string.
        """
        from ..convert import sire_to_rdkit
        from ..legacy.Convert import rdkit_to_smarts
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

        smarts = rdkit_to_smarts(rdkit_mols, map)

        if len(not_water) == 1:
            smarts = [smarts]

        try:
            waters = obj["water"].molecules()
        except Exception:
            from . import SelectorMol

            waters = SelectorMol()

        if len(waters) > 0:
            # need to combine them all together
            s = {}

            for mol, smart in zip(not_water, smarts):
                s[mol.number().value()] = smart

            smarts = []

            if include_hydrogens:
                water = "[#8]1-[H][H]-1"
            else:
                water = "[#8]"

            for mol in obj:
                smarts.append(s.get(mol.number().value(), water))

        if as_search:
            return [f"smarts {x}" for x in smarts]
        else:
            return smarts

    def _selector_to_smiles(
        obj, as_search: bool = True, include_hydrogens: bool = False, map=None
    ):
        """
        Return the molecule views in this container as smiles strings. Include
        hydrogens if 'include_hydrogens' is True. This returns a list
        of smiles strings, in the same order as the views in the container.

        If 'as_search' is True, then this returns the sire search smiles string
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

        if as_search:
            return [f"smiles {x}" for x in smiles]
        else:
            return smiles

    def _selector_view2d(
        obj,
        filename: str = None,
        height: int = 600,
        width: int = 750,
        include_hydrogens: bool = False,
        num_columns: int = 1,
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

        # assign stereochemistry to the rest, and also remove
        # 3D conformers as they mess up the 2D view
        try:
            for r in rdkit_mols:
                # assign the stereochemistry from the structure
                try:
                    from rdkit.Chem.rdmolops import AssignStereochemistryFrom3D

                    AssignStereochemistryFrom3D(r)
                except Exception:
                    # does not matter if this fails
                    pass

                # we also don't want any conformers,
                # as these mess up the 2D view
                r.RemoveAllConformers()
        except Exception:
            try:
                from rdkit.Chem.rdmolops import AssignStereochemistryFrom3D

                AssignStereochemistryFrom3D(rdkit_mols)
            except Exception:
                # does not matter if this fails
                pass

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

    def _to_smarts(
        obj, as_search: bool = False, include_hydrogens: bool = False, map=None
    ):
        """
        Return this molecule view as a smarts string. Include
        hydrogens in 'include_hydrogens' is True.

        If 'as_search' is True, then return this as a sire search
        string
        """
        from ..convert import sire_to_rdkit
        from ..legacy.Convert import rdkit_to_smarts
        from ..base import create_map

        # check if this is a water molecule
        if obj.selected_all():
            try:
                obj = obj["water"]
                if obj.selected_all():
                    if include_hydrogens:
                        smarts = "[#8]1-[H][H]-1"
                    else:
                        smarts = "[#8]"

                    if as_search:
                        return f"smarts {smarts}"
                    else:
                        return smarts
            except Exception:
                pass

        map = create_map(map)

        rdkit_mol = sire_to_rdkit(obj.extract(), map)

        if not include_hydrogens:
            from ..legacy.Convert import rdkit_remove_hydrogens

            rdkit_mol = rdkit_remove_hydrogens(rdkit_mol, map)

        smarts = rdkit_to_smarts(rdkit_mol, map)

        if as_search:
            return f"smarts {smarts}"
        else:
            return smarts

    def _to_smiles(
        obj, as_search: bool = False, include_hydrogens: bool = False, map=None
    ):
        """
        Return this molecule view as a smiles string. Include
        hydrogens in 'include_hydrogens' is True.

        Return this as a sire search string is 'as_search' is True
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
                        smiles = "[H]O[H]"
                    else:
                        smiles = "O"

                    if as_search:
                        return f"smiles {smiles}"
                    else:
                        return smiles
            except Exception:
                pass

        map = create_map(map)

        rdkit_mol = sire_to_rdkit(obj.extract(), map)

        if not include_hydrogens:
            from ..legacy.Convert import rdkit_remove_hydrogens

            rdkit_mol = rdkit_remove_hydrogens(rdkit_mol, map)

        smiles = rdkit_to_smiles(rdkit_mol, map)

        if as_search:
            return f"smiles {smiles}"
        else:
            return smiles

    def _view2d(
        obj,
        filename: str = None,
        height: int = 300,
        width: int = 750,
        include_hydrogens: bool = False,
        num_columns: int = 1,  # included for compatibilty with selector_view2d
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

            # assign the stereochemistry from the structure
            try:
                from rdkit.Chem.rdmolops import AssignStereochemistryFrom3D

                AssignStereochemistryFrom3D(rdkit_mol)
            except Exception:
                # does not matter if this fails
                pass

            # remove conformers, as these mess up the 2D view
            rdkit_mol.RemoveAllConformers()

        if not include_hydrogens:
            from ..legacy.Convert import rdkit_remove_hydrogens

            rdkit_mol = rdkit_remove_hydrogens(rdkit_mol, map)

        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
        from rdkit.Chem.Draw import IPythonConsole

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

    def _selector_to_smarts(obj, *args, **kwargs):
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

    def _selector_to_smarts(obj, *args, **kwargs):
        raise ImportError(
            "You need to install rdkit to be able to generate "
            "smarts strings Do this by typing, e.g. "
            "'mamba install -c conda-forge rdkit' and then restarting "
            "Python and running this script/notebook again."
        )

    def _selector_view2d(obj, *args, **kwargs):
        _view2d(obj, *args, **kwargs)
