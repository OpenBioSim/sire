import pickle

import pytest

import sire as sr


def test_dynamic_attribute_breaks_pickling():
    """
    General, still-open Sire issue: sire_pickle_suite<T> (wrapper/Qt/
    qdatastream.hpp) never overrides pickle_suite::getstate_manages_dict(),
    so any object using it that also has a dynamic Python attribute set on
    it (populating __dict__) fails to pickle with 'Incomplete pickle support
    (__getstate_manages_dict__ not set)', even though the object's own
    QDataStream-based serialisation (sire.stream.save/load) is unaffected.
    """
    from sire.mm import BoreschRestraint, BoreschRestraints

    b = BoreschRestraint(
        receptor=[1574, 1554, 1576],
        ligand=[4, 3, 5],
        r0=sr.u("7.687 A"),
        theta0=[sr.u("1.3031 rad"), sr.u("1.4777 rad")],
        phi0=[sr.u("2.5569 rad"), sr.u("2.9359 rad"), sr.u("1.4147 rad")],
        kr=sr.u("6.2012 kcal mol-1 A-2"),
        ktheta=[sr.u("28.7685 kcal mol-1 rad-2"), sr.u("24.8204 kcal mol-1 rad-2")],
        kphi=[
            sr.u("59.8626 kcal mol-1 rad-2"),
            sr.u("0.7923 kcal mol-1 rad-2"),
            sr.u("55.1775 kcal mol-1 rad-2"),
        ],
    )
    restraints = BoreschRestraints(b)

    # Pickles fine before any dynamic attribute is set.
    pickle.dumps(restraints)

    restraints._some_dynamic_attribute = True

    with pytest.raises(RuntimeError, match="Incomplete pickle support"):
        pickle.dumps(restraints)
