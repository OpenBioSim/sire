from typing import Union as _Union
from typing import List as _List

__all__ = [
    "load",
    "save",
    "save_to_string",
    "expand",
    "tutorial_url",
    "load_test_files",
    "supported_formats",
    "smiles",
    "smarts",
]


class _tutorial_url:
    def __init__(self, value):
        self._value = value
        self.__doc__ = (
            "The base URL for all molecule files used in the tutorial."
        )

    def __str__(self):
        return self._value

    def __repr__(self):
        return self._value

    def startswith(self, value):
        return self._value.startswith(value)


tutorial_url = _tutorial_url("https://sire.openbiosim.org/m")

_range = range


def supported_formats():
    """
    Return a string that describes all of the molecular file formats
    that are supported by Sire
    """
    from .legacy.IO import MoleculeParser

    try:
        return MoleculeParser.supportedFormats()
    except AttributeError:
        return MoleculeParser.supported_formats()


def _create_dir(directory):
    import os

    if not os.path.exists(directory):
        os.makedirs(directory)

    if not os.path.isdir(directory):
        raise IOError(f"{directory} is not a directory!")


def _get_gromacs_dir():
    import os

    if "GROMACS_HOME" in os.environ:
        gromacs_dir = os.environ["GROMACS_HOME"]
        if os.path.exists(gromacs_dir) and os.path.isdir(gromacs_dir):
            return gromacs_dir

    from .config import share_directory

    gromacs_dir = os.path.join(share_directory, "gromacs")

    if os.path.exists(gromacs_dir):
        return gromacs_dir

    # it doesn't exist, so we need to download it
    gromacs_tbz2 = os.path.join(share_directory, "gromacs.tar.bz2")

    if not os.path.exists(gromacs_tbz2):
        try:
            import urllib.request

            urllib.request.urlretrieve(
                f"{tutorial_url}/gromacs.tar.bz2", gromacs_tbz2
            )
        except Exception:
            # we cannot download - just give up
            return None

    if not os.path.exists(gromacs_tbz2):
        return None

    try:
        import tarfile

        t = tarfile.open(gromacs_tbz2, "r|bz2")
        t.extractall(path=share_directory)
    except Exception:
        return None

    if os.path.exists(gromacs_dir):
        return gromacs_dir
    else:
        return None


def _resolve_path(path, directory=".", auto_unzip=True, silent=False):
    import os

    if hasattr(path, "strpath"):
        path = path.strpath

    if hasattr(directory, "strpath"):
        directory = directory.strpath

    if os.path.isdir(path):
        # we need to process this as a trajectory directory - return
        # this as a directory
        return [path]

    elif os.path.exists(path) and os.path.isfile(path):
        if path.endswith(".gz"):
            # unzip the file first
            unzipped = path[0:-3]

            if os.path.exists(unzipped) and os.path.isfile(unzipped):
                if not silent:
                    print(f"Using cached unzipped file '{unzipped}'...")
                return [os.path.abspath(unzipped)]

            _create_dir(directory)
            unzipped = os.path.join(directory, os.path.basename(path)[0:-3])
            if os.path.exists(unzipped) and os.path.isfile(unzipped):
                if not silent:
                    print(f"Using cached unzipped file '{unzipped}'...")
                return [os.path.abspath(unzipped)]

            if auto_unzip:
                if not silent:
                    print(f"Unzipping '{path}'...")

                import gzip
                import shutil

                with gzip.open(path, "rb") as f_in:
                    with open(unzipped, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)

                return [os.path.abspath(unzipped)]
            else:
                return [os.path.abspath(path)]

        elif path.endswith(".bz2"):
            # unzip the file first
            unzipped = path[0:-4]

            if os.path.exists(unzipped) and os.path.isfile(unzipped):
                if not silent:
                    print(f"Using cached unzipped file '{unzipped}'...")
                return [os.path.abspath(unzipped)]

            _create_dir(directory)
            unzipped = os.path.join(directory, os.path.basename(path)[0:-4])
            if os.path.exists(unzipped) and os.path.isfile(unzipped):
                if not silent:
                    print(f"Using cached unzipped file '{unzipped}'...")
                return [os.path.abspath(unzipped)]

            if auto_unzip:
                if not silent:
                    print(f"Unzipping '{path}'...")

                import bz2
                import shutil

                with bz2.open(path, "rb") as f_in:
                    with open(unzipped, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)

                return [os.path.abspath(unzipped)]
            else:
                return [os.path.abspath(path)]

        else:
            return [os.path.abspath(path)]

    if path.startswith("http"):
        # try to download this from the internet
        _create_dir(directory)
        filename = os.path.join(directory, path.split("/")[-1])

        if os.path.exists(filename):
            if os.path.isfile(filename):
                if not silent:
                    print(f"Using cached download of '{path}'...")
                return _resolve_path(
                    filename, directory=directory, silent=silent
                )
            else:
                raise IOError(
                    f"Cannot overwrite {filename} as it is an "
                    "existing directory!"
                )

        if not silent:
            print(f"Downloading from '{path}'...")

        if not filename.endswith(".bz2"):
            # try the bz2 file first
            try:
                import urllib.request

                urllib.request.urlretrieve(f"{path}.bz2", f"{filename}.bz2")
                have_downloaded_file = True
                filename = f"{filename}.bz2"
            except Exception:
                have_downloaded_file = False
        else:
            have_downloaded_file = False

        if not have_downloaded_file:
            try:
                import urllib.request

                urllib.request.urlretrieve(path, filename)
            except Exception as e:
                raise IOError(f"Unable to download '{path}': {e}")

        if os.path.exists(filename) and os.path.isfile(filename):
            return _resolve_path(filename, directory=directory, silent=silent)
        else:
            raise IOError(f"Downloaded file does not exist? {filename}")

    elif len(path) == 4:
        # the first character should be a number
        try:
            int(path[0])
            is_code = True
        except Exception:
            is_code = False

        if is_code:
            code = path.lower()
            # https://files.rcsb.org/download/4hhb.pdb.gz
            return _resolve_path(
                f"https://files.rcsb.org/download/{path}.pdb.gz",
                directory=directory,
                silent=silent,
            )
    elif path.startswith("alphafold:"):
        # alphafold code
        code = path[10:]
        # https://alphafold.ebi.ac.uk/files/AF-" + pdbid + "-F1-model_v1.pdb
        return _resolve_path(
            f"https://alphafold.ebi.ac.uk/files/AF-{code}-F1-model_v3.pdb",
            directory=directory,
            silent=silent,
        )

    # this may be a globbed path
    import glob

    matches = glob.glob(path)

    if len(matches) > 0:
        paths = []
        for match in matches:
            paths += _resolve_path(match, directory=directory, silent=silent)

        return paths

    raise IOError(f"Cannot find file '{path}'")


def expand(base: str, path: _Union[str, _List[str]], *args, **kwargs):
    """Expand the set of paths with the supplied base.

    Args:
     base (str):
         The base to be prepended to all paths

     path (str or list[str]):
         The filename (or names) that will be prepended
         with the base.

     suffix (str):
         A suffix to attach to all files, e.g. ".bz2"

     Returns:
         list[str]:
         The list of expanded filenames or URLs

     Examples:
         >>> expand("https://sire.openbiosim.org/m", "urea.gro", "urea.top")
         ["https://sire.openbiosim.org/m/urea.gro", "https://sire.openbiosim.org/n/urea.top"]

         >>> expand("input", ["ala.top", "ala.crd"])
         ["input/ala.top", "input/ala.crd"]
    """
    if "suffix" in kwargs:
        suffix = kwargs["suffix"]
    else:
        suffix = None

    if type(path) is not list:
        paths = [path]
    else:
        paths = path

    for arg in args:
        paths.append(arg)

    expanded = []

    if base.startswith("http"):

        def join(x, y):
            return f"{x}/{y}"

    else:
        import os

        join = os.path.join

    for path in paths:
        if suffix is None:
            expanded.append(join(base, path))
        else:
            expanded.append(join(base, f"{path}{suffix}"))

    return expanded


def load(
    path: _Union[str, _List[str]],
    *args,
    show_warnings=True,
    silent: bool = False,
    directory: str = ".",
    ignore_topology_frame: bool = False,
    gromacs_path: str = None,
    parallel: bool = True,
    map=None,
    **kwargs,
):
    """
    Load the molecular system at 'path'. This can be a filename
    of a URL. If it is a URL, then the file will be downloaded
    to the current directory and loaded from there.

    Args:
     path (str or list[str]):
        The filename (or names) or the URL or URLS of the molecular
        system to load. This allows multiple paths to be input
        as some molecular file formats split molecular information
        across multiple files. Multiple paths can also be passed
        as multiple arguments to this function.

     show_warnings (bool):
        Whether or not to print out any warnings that are encountered
        when loading your file(s). This is default True, and may lead
        to noisy output. Set `show_warnings=False` to silence this output.

     silent (bool):
        Whether or not to silence all output (including any warnings)

     directory (str):
        Optional directory which will be used when creating any
        files (e.g. as a download from a URL or when unzipping files)

     ignore_topology_frame (bool):
        Ignore any coordinate / frame data coming from the topology file.
        By default, frame data from topology files will be included.
        Setting this to True will ignore that data, meaning that frame
        data will only come from the trajectory files that are loaded.

     gromacs_path (str):
        Path to the directory containing gromacs parameters. If this
        is not set then the gromacs parameters installed with
        sire will be used.

     parallel (bool):
        Whether or not to load files in parallel (using multiple cores).
        You normally do want to do this. Only switch this to False
        if debugging or if you don't want to use all the cores
        in your computer.

    Returns:
        sire.system.System:
        The molecules that have been loaded are returned as
        a sire.system.System

    Examples:
         >>> mols = load("caffeine.pdb")

         >>> mols = load(["ala.crd", "ala.top"])

         >>> mols = load("ala.crd", "ala.top")

         >>> mols = load("https://something")
    """
    if type(path) is not list:
        paths = [path]
    else:
        paths = path

    for arg in args:
        paths.append(arg)

    if silent:
        show_warnings = False

    p = []

    for i in range(0, len(paths)):
        # resolve the paths, downloading as needed
        p += _resolve_path(paths[i], directory=directory, silent=silent)

    paths = p

    if len(paths) == 0:
        raise IOError("No valid files specified. Nothing to load?")

    from .io import load_molecules
    from .base import create_map

    if gromacs_path is None:
        gromacs_path = _get_gromacs_dir()

    m = {
        "GROMACS_PATH": _get_gromacs_dir(),
        "show_warnings": show_warnings,
        "parallel": parallel,
        "ignore_topology_frame": ignore_topology_frame,
    }

    for key in kwargs.keys():
        m[key] = kwargs[key]

    from .base import create_map

    map = create_map(map, m)

    return load_molecules(paths, map=create_map(map))


def _to_legacy_system(molecules):
    """
    Internal function to convert the passed set of molecule views
    into a sire.legacy.System.System
    """
    if hasattr(molecules, "_to_legacy_system"):
        return molecules._to_legacy_system()

    from .legacy.System import System as LegacySystem

    s = LegacySystem()

    if hasattr(molecules, "to_molecule_group"):
        s.add(molecules.to_molecule_group())
    else:
        from .legacy.Mol import MoleculeGroup

        m = MoleculeGroup("all")

        if type(molecules) is list:
            for molecule in molecules:
                m.add(molecule)
        else:
            m.add(molecules)

        s.add(m)

    return s


def save_to_string(
    molecules,
    format: str,
    show_warnings=True,
    silent: bool = False,
    parallel: bool = True,
    map=None,
    **kwargs,
) -> _List[str]:
    """
    Save the passed molecules to an in-memory list of lines.
    This will write the molecule(s) in the format specified
    to memory, thereby avoiding writing any data to a text file

    Note that you must pass in the format, and only a single
    "file" can be written at a time.
    """
    from .base import create_map
    from .legacy.IO import MoleculeParser

    if silent:
        show_warnings = False

    m = {"parallel": parallel, "show_warnings": show_warnings}

    for key in kwargs.keys():
        m[key] = kwargs[key]

    map = create_map(map, m)

    molecules = _to_legacy_system(molecules)

    return MoleculeParser.parse(molecules, format, map=map).lines()


def save(
    molecules,
    filename: str,
    format: _Union[str, _List[str]] = None,
    show_warnings=True,
    silent: bool = False,
    directory: str = ".",
    parallel: bool = True,
    map=None,
    **kwargs,
) -> _List[str]:
    """Save the passed molecules to a file called 'filename'. If the format
    is not specified, then the format will be guessed from the
    filename. If the format is specified, and is a list, then multiple
    files will be written, one for each specified format.

    Args:
     molecules :class:`sire.system.System`,
               :class:`sire.mol.Molecule`,
               List[:class:`sire.mol.Molecule`] etc.)
         The molecule (or molecules) that should be written to the file.
         This can be anything that can be converted to a
         :class:`sire.system.System`, i.e. a single
         :class:`~sire.mol.Molecule` (or :class:`~sire.mol.MoleculeView`),
         or a list of Molecules (or MoleculeViews)

     filename (str):
         The name of the file to which to write the file. Extensions
         will be automatically added if they are needed to match
         the formats of the file (or files) that are written.

     format (str or list(str)):
         The format (or formats) that should be used to write the
         file (or files). If the format isn't specified, then it
         will be guessed from the extension used for `filename`.
         If this doesn't have an extension, then it will be guessed
         based on the formats used to load the molecule originally.
         If it still isn't available, then PDB will be used.

      show_warnings (bool):
         Whether or not to write out any warnings that occur during save

      silent (bool):
         Whether or not to silence all output during the save

      directory (str):
         If supplied, the directory in which to save the files.

      parallel (bool):
         Whether or not to save in parallel (using multiple cores).
         You normally want this switched on, unless you are debugging
         or want to restrict sire to a single core.

    Returns:
         list[str]:
         The absolute paths/name(s) of the files that have been written.

    Examples:
         >>> save(molecules, "molecules.pdb")
         ["/path/to/molecules.pdb"]

         >>> save([mol1, mol2, mol3], "molecules.sdf")
         ["/path/to/molecules.sdf"]

         >>> save(mols, "ala", format=["top", "crd"])
         ["/path/to/ala.top", "/path/to/ala.crd"]
    """
    from .legacy.IO import MoleculeParser
    from .base import create_map

    if silent:
        show_warnings = False

    m = {"parallel": parallel, "show_warnings": show_warnings}

    for key in kwargs.keys():
        m[key] = kwargs[key]

    map = create_map(map, m)

    if hasattr(filename, "strpath"):
        filename = filename.strpath

    if directory is not None:
        if hasattr(directory, "strpath"):
            directory = directory.strpath

        import os

        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = os.path.join(directory, filename)

    if format is not None:
        if type(format) is str:
            format = [format]

        map.set("fileformat", ",".join(format))

    if hasattr(molecules, "_is_trajectory_iterator"):
        # Doing it this way rather that using type(molecules)
        # as type(molecules) randomly fails, and because
        # this way is more pythonic
        if molecules._is_trajectory_iterator():
            # we are saving a trajectory - not just the molecules
            map = molecules._populate_map(map)
            molecules = molecules.current()

    molecules = _to_legacy_system(molecules)

    return MoleculeParser.save(molecules, filename, map=map)


def load_test_files(files: _Union[_List[str], str], *args, map=None):
    """Load the passed files that are part of the unit testing
    and return the resulting molecules. This will cache the files
    into a directory called "../cache" so that downloads can be shared
    between tests. You should only need this function if you
    are writing unit tests.

    Args:
     files (str or list[str])
         The list of files to load from the tutorial website. This
         will automatically add on the tutorial URL and compression suffix

    Returns:
     sire.system.System
         The loaded molecules
    """
    if not type(files) is list:
        files = [files]

    for arg in args:
        files.append(arg)

    import os

    d = os.path.abspath(os.path.curdir)

    if d.endswith("tests"):
        # we are running in the tests directory, so cache downloads here
        cache_dir = os.path.join(d, "cache")
    else:
        d2 = os.path.split(d)[0]
        if d2.endswith("tests"):
            # we are a subdirectory of the parent directory
            cache_dir = os.path.join(d2, "cache")
        else:
            cache_dir = os.path.join(d, "cache")

    files = expand(tutorial_url, files, suffix=".bz2")
    return load(
        files, directory=cache_dir, silent=True, show_warnings=False, map=map
    )


def smiles(
    smiles: str,
    label: str = None,
    labels: str = None,
    smiles_column: str = "smiles",
    labels_column: str = "labels",
    add_hydrogens: bool = True,
    generate_coordinates: bool = True,
    must_sanitize: bool = True,
    map=None,
):
    """
    Return a molecule that has been generated using the passed
    smiles string. This uses rdkit to create the molecule,
    so it must be installed.

    Args:
        smiles: str or list[str] or pandas.Dataframe
            The smiles string to interpret. This can be a single smiles string,
            a list of smiles strings, or a pandas Dataframe containing
            a smiles column and a label column (either called this, or
            use options below to name them yourself)
        label: str
            The label for the molecule being created. This can only
            be a single string. If it is set, then `labels` will be
            ignored.
        labels: str or list[str]
            The label (name) for the molecule that will be created.
            This should be a single string or a list of strings depending
            on 'smiles'. Note that this will be ignored if a
            Dataframe is passed in. Note that if this is not passed in
            then the label will be taken from the smiles string
        smiles_column: str
            The name of the smiles column in the Dataframe (default 'smiles')
        labels_column: str
            The name of the labels column in the Dataframe (default 'labels')
        add_hydrogens: bool (default True)
            Whether or not to automatically add hydrogens. Note that
            not adding hydrogens will automatically disable the
            generation of coordinates.
        generate_coordinates: bool (default True)
            Whether or not to automatically generate 3D coordinates.
            Note that generating coordinates requires that
            hydrogens are automatically added.
        must_sanitize: bool (default True)
            Whether or not all sanity checks must pass when creating
            the molecule. This will ensure that all sanity checks pass,
            and if they don't, then an exception will be raised.
            If this is not True, then sanity checks that failed are
            skipped and silently ignored. It is possible, in this case,
            that a null or malformed molecule may be returned.
        map:
            Property map if you want to put the molecule properties
            into different places

    Returns: sire.mol.Molecule
        The actual molecule
    """
    from .convert import rdkit_to_sire
    from .legacy.Convert import smiles_to_rdkit

    if hasattr(smiles, "to_csv"):
        # convert to a pair of lists from the dataframe
        labels = smiles[[labels_column]]
        smiles = smiles[[smiles_column]]

    elif type(smiles) is not list:
        smiles = [smiles]

    if type(smiles) is list:
        if label is not None:
            labels = label

        if labels is None:
            labels = smiles
        elif type(labels) is not list:
            labels = [labels]

        if len(smiles) != len(labels):
            raise ValueError(
                f"The number of smiles strings {len(smiles)} must match the "
                f"number of labels ({len(labels)})"
            )

    from .base import create_map

    map = create_map(
        map,
        {
            "add_hydrogens": add_hydrogens,
            "generate_coordinates": generate_coordinates,
            "must_sanitize": must_sanitize,
        },
    )

    rdkit_mols = smiles_to_rdkit(smiles, labels, map)

    mols = rdkit_to_sire(rdkit_mols)

    if must_sanitize:
        if len(smiles) == 1:
            mol = mols
            if mol.num_atoms() == 0:
                raise ValueError(
                    "Failed to generate a molecule from the smiles string "
                    f"'{smiles[0]}'. Re-run this function setting "
                    "'must_sanitize' to False if you want to try again, "
                    "ignoring the sanitization steps that failed."
                )
        else:
            empty_mols = []

            for i, mol in enumerate(mols):
                if mol.num_atoms() == 0:
                    empty_mols.append(smiles[i])

            if len(empty_mols) > 0:
                empty_mols = ", ".join(empty_mols)

                raise ValueError(
                    "Failed to generate some molecules from smiles strings. "
                    f"Failed conversions were: [{empty_mols}]. Re-run setting "
                    "'must_sanitize' to False to try to generate the molecule "
                    "ignoring the errors."
                )

    return mols


def smarts(
    smarts: str,
    label: str = None,
    labels: str = None,
    smarts_column: str = "smarts",
    labels_column: str = "labels",
    map=None,
):
    """
    Return a molecule that has been generated using the passed
    smiles string. This uses rdkit to create the molecule,
    so it must be installed.

    Args:
        smiles: str or list[str] or pandas.Dataframe
            The smiles string to interpret. This can be a single smiles string,
            a list of smiles strings, or a pandas Dataframe containing
            a smiles column and a label column (either called this, or
            use options below to name them yourself)
        label: str
            The label for the molecule being created. This can only
            be a single string. If it is set, then `labels` will be
            ignored.
        labels: str or list[str]
            The label (name) for the molecule that will be created.
            This should be a single string or a list of strings depending
            on 'smiles'. Note that this will be ignored if a
            Dataframe is passed in. Note that if this is not passed in
            then the label will be taken from the smiles string
        smiles_column: str
            The name of the smiles column in the Dataframe (default 'smiles')
        labels_column: str
            The name of the labels column in the Dataframe (default 'labels')
        map:
            Property map if you want to put the molecule properties
            into different places

    Returns: sire.mol.Molecule
        The actual molecule
    """
    from .convert import rdkit_to_sire
    from .legacy.Convert import smarts_to_rdkit

    if hasattr(smarts, "to_csv"):
        # convert to a pair of lists from the dataframe
        labels = smiles[[labels_column]]
        smarts = smiles[[smarts_column]]

    elif type(smarts) is not list:
        smarts = [smarts]

    if type(smarts) is list:
        if label is not None:
            labels = label

        if labels is None:
            labels = smarts
        elif type(labels) is not list:
            labels = [labels]

        if len(smarts) != len(labels):
            raise ValueError(
                f"The number of smarts strings {len(smarts)} must match the "
                f"number of labels ({len(labels)})"
            )

    from .base import create_map

    map = create_map(map)

    rdkit_mols = smarts_to_rdkit(smarts, labels, map)

    mols = rdkit_to_sire(rdkit_mols)

    if len(smarts) == 1:
        mol = mols
        if mol.num_atoms() == 0:
            raise ValueError(
                "Failed to generate a molecule from the smarts string "
                f"'{smarts[0]}'."
            )
    else:
        empty_mols = []

        for i, mol in enumerate(mols):
            if mol.num_atoms() == 0:
                empty_mols.append(smarts[i])

        if len(empty_mols) > 0:
            empty_mols = ", ".join(empty_mols)

            raise ValueError(
                "Failed to generate some molecules from smarts strings. "
                f"Failed conversions were: [{empty_mols}]."
            )

    return mols
