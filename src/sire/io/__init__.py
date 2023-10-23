__all__ = ["get_coords_array", "resolve_path"]

from ..legacy import IO as _IO
from .. import use_new_api as _use_new_api

# Imported to ensure that sire.maths.Vector is properly wrapped
from ..maths import Vector as _Vector  # noqa: F401

# Importing the parser sub-module so that it autocompletes when used
from . import parser  # noqa: F401

_use_new_api()


def resolve_path(path, directory=".", auto_unzip=True, silent=False):
    """
    Resolve the passed path into a local path. This will automatically
    download the files associated with a URL path to local files
    that are placed into the specified `directory`.

    This returns the full path to the locally resolved file.

    Args:
        path: The full path that needs resolving. If this is a URL
              then it will be downloaded (either the original URL,
              of a bzip2 version if the original doesn't exist).
              If this is a PDB or alphafold code, then the corresponding
              PDB or alphafold structure will be downloaded.

        directory: The full path to the directory you want any downloaded
                   files to be downloaded into.

        auto_unzip: Automatically unzip any compressed files after download
                    (or if they are local). Note that the original
                    zipped file is not removed.

        silent: Silence any output printed to the screen when resolving
                the path.

    Returns:
        The list of filenames of all the files that are downloaded / resolved
    """
    from .._load import _resolve_path

    return _resolve_path(
        path=path, directory=directory, auto_unzip=auto_unzip, silent=silent
    )


def get_coords_array(mol, units=None, map=None):
    """
    Return the coordinates of the passed molecule view as a
    numpy array of shape (natoms,3). Specify the length
    units to use, and optionally pass in a map to find
    the coordinates property
    """
    import numpy as np

    if units is None:
        from ..units import angstrom

        units = angstrom

    from ..base import create_map

    map = create_map(map)

    if hasattr(mol, "to_molecule_group"):
        mol = mol.to_molecule_group()

    coords = _IO.getCoordsArray(mol, units, map)

    natoms = int(len(coords) / 3)

    return np.reshape(np.asarray(coords, dtype=float), (natoms, 3))


def load_molecules(*args, map=None, **kwargs):
    from ..legacy.IO import load_molecules as _load_molecules
    from ..system import System
    from ..base import create_map

    map = create_map(map)

    mols = System(_load_molecules(*args, map=map, **kwargs))

    try:
        mols.add_shared_property("space", mols.property("space"))
    except Exception:
        from ..vol import Cartesian

        mols.add_shared_property("space", Cartesian())

    try:
        mols.add_shared_property("time", mols.property("time"))
    except Exception:
        from ..units import picosecond
        from ..base import wrap

        mols.add_shared_property("time", wrap(0 * picosecond))

    if map.specified("make_whole"):
        if map["make_whole"]:
            mols.make_whole()

    return mols
