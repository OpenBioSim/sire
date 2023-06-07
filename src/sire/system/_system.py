__all__ = ["System"]


class System:
    """
    This class holds a complete system of molecules. This is normally
    created by reading in molecules from an input file.

    This provides a "MoleculeView"-style interface to the molecules,
    acting very similarly to a sire.mol.SelectorMol object.

    You can convert this to a sire.mol.SelectorMol object by
    calling the System.molecules() function.
    """

    def __init__(self, system=None):
        from ..legacy.System import System as _System

        if system is None:
            self._system = _System()
        else:
            if _System not in type(system).mro():
                raise TypeError(
                    "You can only construct from a sire.legacy.System.System, "
                    f"not a {type(system)}"
                )

            if type(system) == System:
                self._system == system._system
            else:
                self._system = system

        self._molecules = None

    @staticmethod
    def is_system(obj):
        """
        Return whether the passed object is a System class
        (either a new or legacy System)
        """
        from ..legacy.System import System as _System

        return type(obj) == System or type(obj) == _System

    def __copy__(self):
        other = System()
        other._system = self._system.clone()
        other._molecules = None
        return other

    def __deepcopy__(self, memo):
        return self.__copy__()

    def __str__(self):
        return str(self._system)

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, key):
        return self.molecules()[key]

    def __iadd__(self, molecules):
        self.add(molecules)
        return self

    def __add__(self, molecules):
        ret = self.__copy__()
        ret.add(molecules)
        return ret

    def __radd__(self, molecules):
        return self.__add__(molecules)

    def __isub__(self, molecules):
        self.remove(molecules)
        return self

    def __sub__(self, molecules):
        ret = self.__copy__()
        ret.remove(molecules)
        return ret

    def clone(self):
        """Return a copy (clone) of this System"""
        s = System()
        s._system = self._system.clone()
        s._molecules = None
        return s

    def count(self):
        """Return the number of items in this System"""
        return self.__len__()

    def size(self):
        """Return the number of items in this System"""
        return self.__len__()

    def __len__(self):
        return len(self.molecules())

    def num_atoms(self):
        """Return the number of atoms in this System"""
        return self._system.num_atoms()

    def num_residues(self):
        """Return the number of residues in this System"""
        return self._system.num_residues()

    def num_chains(self):
        """Return the number of chains in this System"""
        return self._system.num_chains()

    def num_segments(self):
        """Return the number of segments in this System"""
        return self._system.num_segments()

    def num_molecules(self):
        """Return the number of molecules in this System"""
        return self._system.num_molecules()

    def find(self, views):
        """Return the index(es) of the molecule(s) that are in `views`"""
        return self.molecules().find(views)

    def names(self):
        """Return the names of all of the molecules in this System"""
        return self.molecules().names()

    def numbers(self):
        """Return the numbers of all of the molecules in this System"""
        return self.molecules().numbers()

    def num_frames(self):
        """Return the number of trajectory frames for this System"""
        return self._system.num_frames()

    def load_frame(self, i):
        """Load the ith frame into this System"""
        self._system.load_frame(i)
        self._molecules = None

    def save_frame(self, i=None, map=None):
        """Save the current coordinates to the ith frame of this System.
        If i is not specfied then this adds the frame onto the
        end of the trajectory
        """
        from ..base import create_map

        map = create_map(map)

        if i is None:
            self._system.save_frame(map=map)
        else:
            self._system.save_frame(i, map=map)

        self._molecules = None

    def delete_frame(self, i):
        """Delete the ith frame from the trajectory"""
        self._system.delete_frame(i)
        self._molecules = None

    def to_molecule_group(self):
        """Return this System converted to a sire.mol.MoleculeGroup"""
        return self.molecules().to_molecule_group()

    def molecules(self, *args, **kwargs):
        """Return this System converted to a sire.mol.SelectorMol.
        You can pass in arguments to search or index so that you
        limit the number of molecules returned.
        """
        if self._molecules is not None:
            return self._molecules.molecules(*args, **kwargs)

        import sire.mol

        self._molecules = sire.mol.SelectorMol(self._system)

        if self._molecules.num_atoms() != self._system.num_atoms():
            # oh dear - this is an edge case where the System does
            # not contain complete molecules. We need to extract
            # the molecules and re-add them
            raise NotImplementedError(
                "sire.system.System does not yet support Systems that hold "
                "partial molecules!. Let us know that you have hit this "
                "bug and we will add support."
            )

        return self.molecules(*args, **kwargs)

    def segments(self, *args, **kwargs):
        """Return all segments in this System (or those that match
        the passed index, if supplied)
        """
        return self.molecules().segments(*args, **kwargs)

    def chains(self, *args, **kwargs):
        """Return all chains in this System (or those that match
        the passed index, if supplied)
        """
        return self.molecules().chains(*args, **kwargs)

    def residues(self, *args, **kwargs):
        """Return all residues in this System (or those that match
        the passed index, if supplied)
        """
        return self.molecules().residues(*args, **kwargs)

    def atoms(self, *args, **kwargs):
        """Return all atoms in this System (or those that match
        the passed index, if supplied)
        """
        return self.molecules().atoms(*args, **kwargs)

    def bonds(self, *args, **kwargs):
        """Return all bonds in this System (or those that match
        the passed index, if supplied)
        """
        return self.molecules().bonds(*args, **kwargs)

    def angles(self, *args, **kwargs):
        """Return all angles in this System (or those that match
        the passed index, if supplied)
        """
        return self.molecules().angles(*args, **kwargs)

    def dihedrals(self, *args, **kwargs):
        """Return all dihedrals in this System (or those that match
        the passed index, if supplied)
        """
        return self.molecules().dihedrals(*args, **kwargs)

    def impropers(self, *args, **kwargs):
        """Return all impropers in this System (or those that match
        the passed index, if supplied)
        """
        return self.molecules().impropers(*args, **kwargs)

    def molecule(self, *args, **kwargs):
        """Return the molecule that matches the passed index/search"""
        return self.molecules().molecule(*args, **kwargs)

    def segment(self, *args, **kwargs):
        """Return the segment that matches the passed index/search"""
        return self.molecules().segment(*args, **kwargs)

    def chain(self, *args, **kwargs):
        """Return the chain that matches the passed index/search"""
        return self.molecules().chain(*args, **kwargs)

    def residue(self, *args, **kwargs):
        """Return the residue that matches the passed index/search"""
        return self.molecules().residue(*args, **kwargs)

    def atom(self, *args, **kwargs):
        """Return the atom that matches the passed index/search"""
        return self.molecules().atom(*args, **kwargs)

    def bond(self, *args, **kwargs):
        """Return the bond that matches the passed index/search"""
        return self.molecules().bond(*args, **kwargs)

    def angle(self, *args, **kwargs):
        """Return the angle that matches the passed index/search"""
        return self.molecules().angle(*args, **kwargs)

    def dihedral(self, *args, **kwargs):
        """Return the dihedral that matches the passed index/search"""
        return self.molecules().dihedral(*args, **kwargs)

    def improper(self, *args, **kwargs):
        """Return the improper that matches the passed index/search"""
        return self.molecules().improper(*args, **kwargs)

    def trajectory(self, *args, **kwargs):
        """Return an iterator over the trajectory of frames for this System"""
        return self.molecules().trajectory(*args, **kwargs)

    def minimisation(self, map=None):
        """
        Return a Minimisation object that can be used to minimise the energy
        of the molecule(s) in this view.
        """
        from ..mol import Minimisation

        return Minimisation(self, map=map)

    def dynamics(self, *args, **kwargs):
        """
        Return a Dynamics object that can be used to perform
        dynamics of the molecule(s) in this view
        """
        from ..mol import _dynamics

        return _dynamics(self, *args, **kwargs)

    def energy(self, *args, **kwargs):
        """Calculate and return the energy of this System
        (or of the matching index/search subset of this System)
        """
        return self.molecules().energy(*args, **kwargs)

    def energies(self, *args, **kwargs):
        """Calculate and return the individual energies of the
        contents of this System (or of the matching index/search
        subset of this System)
        """
        return self.molecules().energies(*args, **kwargs)

    def charge(self, *args, **kwargs):
        """Return the total charge of this System (or of the matching
        index/search subset of this System)
        """
        return self.molecules().charge(*args, **kwargs)

    def mass(self, *args, **kwargs):
        """Return the total mass of this System (or of the matching
        index/search subset of this System)
        """
        return self.molecules().mass(*args, **kwargs)

    def coordinates(self, *args, **kwargs):
        """Return the center of geometry of this System (or of the matching
        index/search subset of this System)
        """
        return self.molecules().coordinates(*args, **kwargs)

    def space(self, map=None):
        """
        Return the space used for this system
        """
        try:
            if map is None:
                return self._system.property("space")
            else:
                from ..base import create_map

                map = create_map(map)
                return self._system.property(map["space"])
        except Exception:
            from ..vol import Cartesian

            return Cartesian()

    def time(self, map=None):
        """
        Return the current system time
        """
        try:
            if map is None:
                return self._system.property("time")
            else:
                from ..base import create_map

                map = create_map(map)
                return self._system.property(map["time"])
        except Exception:
            from ..units import picosecond

            return 0 * picosecond

    def set_space(self, space, map=None):
        """
        Set the space to be used to hold all of the molecules
        in this system
        """
        from ..legacy.Vol import Space

        if not issubclass(type(space), Space):
            raise TypeError(
                "You can only set the space to a type derived from sire.vol.Space, "
                "e.g. sire.vol.PeriodicBox, sire.vol.TriclinicBox or "
                f"sire.vol.Cartesian. You cannot use a {type(space)}."
            )

        if map is None:
            self._system.set_property("space", space)
        else:
            from ..base import create_map

            map = create_map(map)

            space_property = map["space"]

            if space_property.has_source():
                self._system.set_property(space_property.source(), space)

    def set_time(self, time, map=None):
        """
        Set the current time for the system
        """
        from ..units import picosecond
        from ..base import wrap

        if time == 0:
            time = wrap(0 * picosecond)
        else:
            if not hasattr(time, "has_same_units"):
                raise TypeError(
                    "You can only set the time to a value with units 'time', e.g. "
                    f"5 * sire.units.picosecond. YOu cannot use a {type(time)}."
                )

            if not time.has_same_units(picosecond):
                raise TypeError(
                    "You can only set the time to a value of units time. You "
                    f"cannot use {time}."
                )

            time = wrap(time)

        if map is None:
            self._system.set_property("time", time)
        else:
            from ..base import create_map

            map = create_map(map)

            time_property = map["time"]

            if time_property.has_source():
                self._system.set_property(time_property.source(), time)

    def evaluate(self, *args, **kwargs):
        """Return an evaluator for this Systme (or of the matching
        index/search subset of this System)"""
        return self.molecules().evaluate(*args, **kwargs)

    def property(self, *args, **kwargs):
        """Return the System property that matches the passed key/arguments"""
        return self._system.property(*args, **kwargs)

    def set_property(self, *args, **kwargs):
        """Set the System property according to the passed arguments"""
        self._system.set_property(*args, **kwargs)

    def properties(self):
        """Return all of the System-level properties of this System"""
        return self._system.properties()

    def property_keys(self):
        """Return the keys (IDs) of all of the System-level properties
        of this System
        """
        return self._system.property_keys()

    def cursor(self):
        """
        Return a sire.mol.Cursor that can be used to edit
        the molecules in this System
        """
        from ..mol._cursor import CursorsM

        return CursorsM(self)

    def smiles(self, *args, **kwargs):
        """
        Return the molecule views in this container as smiles strings. Include
        hydrogens in 'include_hydrogens' is True. This returns a list
        of smiles strings, in the same order as the views in the container
        """
        return self.molecules().smiles(*args, **kwargs)

    def view(self, *args, **kwargs):
        """
        View this System (or the matching index/search subset)
        via a nglview viewer. Only works in an interactive
        notebook session, e.g. in a Jupyter notebook
        """
        return self.molecules().view(*args, **kwargs)

    def view2d(self, *args, **kwargs):
        """
        Create a 2D representation of the molecules in this system.
        If 'filename' is set then this will be written to that file. Otherwise
        this will be returned for visualisation in a jupyter notebook.
        """
        return self.molecules().view2d(*args, **kwargs)

    def add(self, molecules):
        """
        Add the passed molecules to this system. This will only
        add molecules that don't already exist in this system.
        """
        if type(molecules) is list:
            for molecule in molecules:
                self.add(molecule)
            return

        from ..legacy.Mol import MGName

        mgid = MGName("all")

        if hasattr(molecules, "molecules"):
            # rather convoluted way to get a 'Molecules' object...
            mols = molecules.molecules().to_molecule_group().molecules()
            self._system.add_if_unique(mols, mgid)
        else:
            self._system.add_if_unique(molecules, mgid)

        self._molecules = None

    def remove(self, molecules):
        """
        Remove the passed molecules from this system.
        """
        if type(molecules) is list:
            for molecule in molecules:
                self.remove(molecule)
            return

        molnums = None

        if hasattr(molecules, "molecules"):
            molnums = molecules.molecules().mol_nums()
        else:
            molnums = [molecule.molecule().number()]

        for molnum in molnums:
            self._system.remove(molnum)

        self._molecules = None

    def update(self, molecules):
        """
        Update the molecules in this system so that they have
        the same versions and data as the new molecules contained
        in 'value'
        """
        if type(molecules) is list:
            for molecule in molecules:
                self.update(molecule)
            return

        if hasattr(molecules, "molecules"):
            self._system.update(molecules.molecules())
        else:
            self._system.update(molecules)

        self._molecules = None

    def apply(self, *args, **kwargs):
        """
        Apply (perform / map) the passed function (with associated arguments)
        on all of the molecules in this System
        """
        return self.molecules().apply(*args, **kwargs)

    def apply_reduce(self, *args, **kwargs):
        """
        Apply (perform / map) the passed function (with associated arguments)
        on all of the molecules in this System, reducing the result
        using the passed reduction function
        """
        return self.molecules().apply_reduce(*args, **kwargs)
