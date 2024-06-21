__all__ = ["create_engine", "emle", "zero_charge"]

from .. import use_new_api as _use_new_api

_use_new_api()

from ..legacy import Convert as _Convert

from ._emle import emle
from ._utils import _zero_charge as zero_charge


def create_engine(
    mols,
    qm_atoms,
    py_object,
    callback=None,
    cutoff="7.5A",
    neighbourlist_update_frequency=20,
    redistribute_charge=False,
    map=None,
):
    """
    Create a QM engine to that can be used for QM/MM simulations with sire.mol.dynamics.

    Parameters
    ----------

    mols : sire.system.System
        The molecular system.

    qm_atoms : str, int, list, molecule view/collection etc.
        Any valid search string, atom index, list of atom indicies,
        or molecule view/container that can be used to select
        qm_atoms from 'mols'.

    py_object : object
        The Python object that will contains the callback for the QM calculation.
        This can be a class instance with a "callback" method, or a callable.

    callback : str, optional, default=None
        The name of the callback. If None, then the py_object is assumed to
        be a callable, i.e. it is itself the callback.

    cutoff : str or sire.legacy.Units.GeneralUnit, optional, default="7.5A"
        The cutoff to use for the QM/MM calculation.

    neighbourlist_update_frequency : int, optional, default=20
        The frequency with which to update the neighbourlist.

    redistribute_charge : bool
        Whether to redistribute charge of the QM atoms to ensure that the total
        charge of the QM region is an integer. Excess charge is redistributed
        over the non QM atoms within the residues involved in the QM region.

    Returns
    -------

    engine : sire.legacy.Convert._SireOpenMM.PyQMEngine
        The QM engine.
    """

    from ..base import create_map as _create_map
    from ..mol import selection_to_atoms as _selection_to_atoms
    from ..system import System as _System
    from ..legacy import Units as _Units
    from ..units import angstrom as _angstrom
    from .. import u as _u

    if not isinstance(mols, _System):
        raise TypeError("mols must be a of type 'sire.System'")

    # Clone the system.
    mols = mols.clone()

    try:
        qm_atoms = _selection_to_atoms(mols, qm_atoms)
    except:
        raise ValueError("Unable to select 'qm_atoms' from 'mols'")

    if callback is not None:
        if not isinstance(callback, str):
            raise TypeError("'callback' must be of type 'str'")
        if not hasattr(py_object, callback):
            raise ValueError(f"'py_object' does not have a method called '{callback}'.")
    else:
        callback = ""

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

    if not isinstance(redistribute_charge, bool):
        raise TypeError("'redistribute_charge' must be of type 'bool'")

    if map is not None:
        if not isinstance(map, dict):
            raise TypeError("'map' must be of type 'dict'")
    map = _create_map(map)

    # Create the EMLE engine.
    engine = _Convert.PyQMEngine(
        py_object,
        callback,
        cutoff,
        neighbourlist_update_frequency,
    )

    from ._utils import (
        _check_charge,
        _create_qm_mol_to_atoms,
        _configure_engine,
        _create_merged_mols,
        _get_link_atoms,
    )

    # Check that the charge of the QM region is integer valued.
    _check_charge(mols, qm_atoms, map, redistribute_charge)

    # Get the mapping between molecule numbers and QM atoms.
    qm_mol_to_atoms = _create_qm_mol_to_atoms(qm_atoms)

    # Get link atom information.
    mm1_to_qm, mm1_to_mm2, bond_scale_factors, mm1_indices = _get_link_atoms(
        mols, qm_mol_to_atoms, map
    )

    # Configure the engine.
    engine = _configure_engine(
        engine, mols, qm_atoms, mm1_to_qm, mm1_to_mm2, bond_scale_factors, map
    )

    # Create the merged molecule.
    qm_mols = _create_merged_mols(qm_mol_to_atoms, mm1_indices, map)

    # Update the molecule in the system.
    mols.update(qm_mols)

    # Bind the system as a private attribute of the engine.
    engine._mols = mols

    return mols, engine
