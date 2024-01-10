__all__ = ["emle"]

from ..legacy import Convert as _Convert

from .. import use_new_api as _use_new_api

_use_new_api()

_EMLEEngine = _Convert._SireOpenMM.EMLEEngine

try:
    from emle.emle import EMLECalculator as _EMLECalculator

    _has_emle = True
except:
    _has_emle = False


def emle(
    mols,
    calculator,
    qm_index,
    cutoff="7.5A",
    neighbourlist_update_frequency=20,
    map=None,
):
    """
    Create an EMLEEngine object to allow QM/MM simulations using sire.mol.dynamics.

    Parameters
    ----------

    mols : sire.system.System
        The molecular system.

    calculator : emle.emle.EMLECalculator
        The EMLECalculator object to use for elecotrostatic embedding calculations.

    qm_index : int
        The index of the QM molecule in the system.

    cutoff : str or sire.legacy.Units.GeneralUnit, optional, default="7.5A"
        The cutoff to use for the QM/MM calculation.

    neighbourlist_update_frequency : int, optional, default=20
        The frequency with which to update the neighbourlist.

    Returns
    -------

    engine : sire.legacy.Convert._SireOpenMM.EMLEEngine
        The EMLEEngine object.
    """
    if not _has_emle:
        raise ImportError(
            "Could not import emle. Please install emle-egine and try again."
        )

    from ..base import create_map as _create_map
    from ..system import System as _System
    from ..legacy import Units as _Units
    from ..units import angstrom as _angstrom
    from .. import u as _u

    if not isinstance(mols, _System):
        raise TypeError("mols must be a of type 'sire.System'")

    if not isinstance(calculator, _EMLECalculator):
        raise TypeError("'calculator' must be a of type 'emle.emle.EMLECalculator'")

    if not isinstance(qm_index, int):
        raise TypeError("'qm_index' must be of type 'int'")

    try:
        qm_mol = mols[qm_index]
    except:
        raise ValueError(f"qm_index must be in range [0, {len(mols)})")

    if not isinstance(cutoff, (str, _Units.GeneralUnit)):
        raise TypeError(
            "cutoff must be of type 'str' or 'sire.legacy.Units.GeneralUnit'"
        )

    if isinstance(cutoff, str):
        try:
            cutoff = _u(cutoff)
        except:
            raise ValueError("Unable to parse cutoff as a GeneralUnit")

    if not cutoff.has_same_units(_angstrom):
        raise ValueError("'cutoff' must be in units of length")

    if not isinstance(neighbourlist_update_frequency, int):
        raise TypeError("'neighbourlist_update_frequency' must be of type 'int'")

    if neighbourlist_update_frequency < 0:
        raise ValueError("'neighbourlist_update_frequency' must be >= 0")

    if map is not None:
        if not isinstance(map, dict):
            raise TypeError("'map' must be of type 'dict'")
    map = _create_map(map)

    # Create the EMLEEngine.
    engine = _EMLEEngine(
        calculator,
        cutoff,
        neighbourlist_update_frequency,
    )

    from ._utils import _configure_engine

    # Configure the engine.
    engine = _configure_engine(engine, mols, qm_mol, map)

    # Get the QM property flag.
    qm_propname = map["is_qm"]

    # Check for existing QM molecules. Currently we only support a single
    # molecule.
    try:
        qm_mols = mols[f"property {qm_propname}"].molecules()
    except:
        qm_mols = []

    if len(qm_mols) > 0:
        raise ValueError("This system already contains a QM molecule!")

    # Create a cursor.
    c = qm_mol.cursor()

    # Flag the molecule as QM.
    c[qm_propname] = True

    # Commit the changes and update.
    qm_mol = c.commit()
    mols.update(qm_mol)

    return mols, engine
