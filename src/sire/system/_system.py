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
            from ..vol import Cartesian
            from ..units import picosecond
            from ..base import wrap
            from ..legacy.Mol import MoleculeGroup

            self._system = _System()
            self._system.add_shared_property("space", Cartesian())
            self._system.add_shared_property("time", wrap(0 * picosecond))
            self._system.add(MoleculeGroup("all"))
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

    def _to_legacy_system(self):
        """
        Internal function used to convert this back to a legacy system
        """
        return self._system

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

    def make_whole(self, map=None):
        """
        Make all of the molecules in this system whole. This
        maps each molecule into the current space, such that no
        molecule is broken across a periodic box boundary
        """
        if map is None:
            self._system.make_whole()
        else:
            from ..base import create_map

            self._system.make_whole(map=create_map(map))

        self._molecules = None

    def num_frames(self, map=None):
        """Return the number of trajectory frames for this System"""
        from ..base import create_map

        return self._system.num_frames(map=create_map(map))

    def load_frame(self, i, map=None):
        """Load the ith frame into this System"""
        from ..base import create_map

        self._system.load_frame(i, map=create_map(map))
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

    def delete_frame(self, i, map=None):
        """Delete the ith frame from the trajectory"""
        from ..base import create_map

        self._system.delete_frame(i, map=create_map(map))
        self._molecules = None

    def delete_all_frames(self, map=None):
        """Delete all the frames from the trajectory"""
        from ..base import create_map

        self._system.delete_all_frames(map=create_map(map))
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
        """
        Return an iterator over the trajectory of frames of this view.

        align:
          Pass in a selection string to select atoms against which
          every frame will be aligned. These atoms will be moved
          to the center of the periodic box (if a periodic box
          is used). If "True" is passed, then this will attempt
          to align *ALL* of the coordinates in the view.

          You can also choose to pass in a molecular container,
          and it will align against the atoms in that container,
          assuming they are contained in this view. If not, then
          you need to supply a mapping that maps from the
          atoms in the align container, to the atoms in this view.

        frame:
          The frame of the trajectory against which the alignment
          should be based. For example, `frame=3` would align based
          on the coordinates of the aligned atoms in frame 3 of
          the trajectory. If this is `None` (the default) then the
          first frame will be used.

        mapping: AtomMapping
            An AtomMapping object that maps from atoms in the alignment
            container to atoms in this view. You only need to supply
            this if all of the alignment atoms are not contained
            in this view.

        smooth:
            Pass in the number of frames to smooth (average) the view
            over. If 'True' is passed, then the recommended number
            of frames will be averaged over

        wrap: bool
            Whether or not to wrap the coordinates into the periodic box
        """
        from ..mol._trajectory import TrajectoryIterator

        return TrajectoryIterator(self, *args, **kwargs)

    def minimisation(self, *args, **kwargs):
        """
        Return a Minimisation object that can be used to perform
        minimisation of the molecule(s) in this System

        cutoff: Length
            The size of the non-bonded cutoff

        cutoff_type: str
            The type of cutoff to use, e.g. "PME", "RF" etc.
            See https://sire.openbiosim.org/cheatsheet/openmm.html#choosing-options
            for the full list of options

        constraint: str
            The type of constraint to use for bonds and/or angles, e.g.
            `h-bonds`, `bonds` etc.
            See https://sire.openbiosim.org/cheatsheet/openmm.html#choosing-options
            for the full list of options. This will be automatically
            guessed from the timestep if it isn't set.

        schedule: sire.cas.LambdaSchedule
            The schedule used to control how perturbable forcefield parameters
            should be morphed as a function of lambda. If this is not set
            then a sire.cas.LambdaSchedule.standard_morph() is used.

        lambda_value: float
            The value of lambda at which to run minimisation. This only impacts
            perturbable molecules, whose forcefield parameters will be
            scaled according to the lambda schedule for the specified
            value of lambda.

        swap_end_states: bool
            Whether or not to swap the end states. If this is True, then
            the perturbation will run from the perturbed back to the
            reference molecule (the perturbed molecule will be at lambda=0,
            while the reference molecule will be at lambda=1). This will
            use the coordinates of the perturbed molecule as the
            starting point.

        ignore_perturbations: bool
            Whether or not to ignore perturbations. If this is True, then
            the perturbation will be ignored, and the simulation will
            be run using the properties of the reference molecule
            (or the perturbed molecule if swap_end_states is True). This
            is useful if you just want to run standard molecular dynamics
            of the reference or perturbed states.

        shift_delta: length
            The shift_delta parameter that controls the electrostatic
            and van der Waals softening potential that smooths the
            creation and deletion of ghost atoms during a potential.
            This defaults to 2.0 A.

        coulomb_power: int
            The coulomb power parmeter that controls the electrostatic
            softening potential that smooths the creation and deletion
            of ghost atoms during a potential. This defaults to 0.

        restraints: sire.mm.Restraints or list[sire.mm.Restraints]
            A single set of restraints, or a list of sets of
            restraints that will be applied to the atoms during
            the simulation.

        fixed: molecule(s) view, search string, int, list[int] etc
            Anything that can be used to identify the atom or atoms
            that should be fixed in place during the simulation. These
            atoms will not be moved by minimisation.

        device: str or int
            The ID of the GPU (or accelerator) used to accelerate
            minimisation. This would be CUDA_DEVICE_ID or similar
            if CUDA was used. This can be any valid OpenMM device string

        map: dict
            A dictionary of additional options. Note that any options
            set in this dictionary that are also specified via one of
            the arguments above will be overridden by the argument
            value
        """
        from ..mol import _minimisation

        return _minimisation(self, *args, **kwargs)

    def dynamics(self, *args, **kwargs):
        """
        Return a Dynamics object that can be used to perform
        dynamics on the molecule(s) in this System

        cutoff: Length
            The size of the non-bonded cutoff

        cutoff_type: str
            The type of cutoff to use, e.g. "PME", "RF" etc.
            See https://sire.openbiosim.org/cheatsheet/openmm.html#choosing-options
            for the full list of options

        timestep: time
            The size of the dynamics timestep

        save_frequency: time
            The amount of simulation time between saving energies and frames.
            This can be overridden using `energy_frequency` or `frame_frequency`,
            or by these options in individual dynamics runs. Set this
            to zero if you don't want any saves.

        energy_frequency: time
            The amount of time between saving energies. This overrides the
            value in `save_frequency`. Set this to zero if you don't want
            to save energies during the trajectory. This can be overridden
            by setting energy_frequency during an individual run.

        frame_frequency: time
            The amount of time between saving frames. This overrides the
            value in `save_frequency`. Set this to zero if you don't want
            to save frames during the trajectory. This can be overridden
            by setting frame_frequency during an individual run.

        save_velocities: bool
            Whether or not to save velocities when saving trajectory frames
            during the simulation. This defaults to False, as velocity
            trajectories aren't often needed, and they double the amount
            of space that is required for a trajectory.

        constraint: str
            The type of constraint to use for bonds and/or angles, e.g.
            `h-bonds`, `bonds` etc.
            See https://sire.openbiosim.org/cheatsheet/openmm.html#choosing-options
            for the full list of options. This will be automatically
            guessed from the timestep if it isn't set.

        schedule: sire.cas.LambdaSchedule
            The schedule used to control how perturbable forcefield parameters
            should be morphed as a function of lambda. If this is not set
            then a sire.cas.LambdaSchedule.standard_morph() is used.

        lambda_value: float
            The value of lambda at which to run dynamics. This only impacts
            perturbable molecules, whose forcefield parameters will be
            scaled according to the lambda schedule for the specified
            value of lambda.

        swap_end_states: bool
            Whether or not to swap the end states. If this is True, then
            the perturbation will run from the perturbed back to the
            reference molecule (the perturbed molecule will be at lambda=0,
            while the reference molecule will be at lambda=1). This will
            use the coordinates of the perturbed molecule as the
            starting point.

        ignore_perturbations: bool
            Whether or not to ignore perturbations. If this is True, then
            the perturbation will be ignored, and the simulation will
            be run using the properties of the reference molecule
            (or the perturbed molecule if swap_end_states is True). This
            is useful if you just want to run standard molecular dynamics
            of the reference or perturbed states.

        temperature: temperature
            The temperature at which to run the simulation. A
            microcanonical (NVE) simulation will be run if you don't
            specify the temperature.

        pressure: pressure
            The pressure at which to run the simulation. A
            microcanonical (NVE) or canonical (NVT) simulation will be
            run if the pressure is not set.

        shift_delta: length
            The shift_delta parameter that controls the electrostatic
            and van der Waals softening potential that smooths the
            creation and deletion of ghost atoms during a potential.
            This defaults to 2.0 A.

        coulomb_power: int
            The coulomb power parmeter that controls the electrostatic
            softening potential that smooths the creation and deletion
            of ghost atoms during a potential. This defaults to 0.

        restraints: sire.mm.Restraints or list[sire.mm.Restraints]
            A single set of restraints, or a list of sets of
            restraints that will be applied to the atoms during
            the simulation.

        fixed: molecule(s) view, search string, int, list[int] etc
            Anything that can be used to identify the atom or atoms
            that should be fixed in place during the simulation. These
            atoms will not be moved by dynamics.

        device: str or int
            The ID of the GPU (or accelerator) used to accelerate
            the simulation. This would be CUDA_DEVICE_ID or similar
            if CUDA was used. This can be any valid OpenMM device string

        precision: str
            The desired precision for the simulation (e.g. `single`,
            `mixed` or `double`)

        map: dict
            A dictionary of additional options. Note that any options
            set in this dictionary that are also specified via one of
            the arguments above will be overridden by the argument
            value
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

        self._molecules = None

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

        self._molecules = None

    def energy_trajectory(
        self, to_pandas: bool = True, to_alchemlyb: bool = False, map=None
    ):
        """
        Return the energy trajectory for this System. This is the history
        of energies evaluate during any dynamics runs. It could include
        energies calculated at different values of lambda.

        Parameters
        ----------

        to_pandas: bool
            Whether or not to return the energy trajectory as a pandas DataFrame.

        to_alchemlyb: bool
            Whether or not to return the energy trajectory as a pandas DataFrame
            that is formatted to usable in alchemlyb
        """
        from ..base import create_map

        map = create_map(map)

        traj_propname = map["energy_trajectory"]

        try:
            traj = self._system.property(traj_propname)
        except Exception:
            traj = None

        if traj is not None:
            if traj.what() != "SireMaths::EnergyTrajectory":
                if traj_propname.has_value():
                    raise TypeError(
                        f"You cannot force the use of a {type(traj)} "
                        "as an EnergyTrajectory"
                    )

                traj = None

        if traj is None:
            # we need to create this trajectory
            from ..maths import EnergyTrajectory

            self._system.set_property(
                traj_propname.source(), EnergyTrajectory()
            )

            traj = self._system.property(traj_propname)

        if to_pandas or to_alchemlyb:
            return traj.to_pandas(to_alchemlyb=to_alchemlyb)
        else:
            return traj

    def set_energy_trajectory(self, trajectory, map=None):
        """
        Set the energy trajectory to the passed value
        """
        from ..base import create_map

        map = create_map(map)

        traj_propname = map["energy_trajectory"]

        if traj_propname.has_value():
            return

        if trajectory.what() != "SireMaths::EnergyTrajectory":
            raise TypeError(
                f"You cannot set a {type(trajectory)} as an "
                "energy trajectory!"
            )

        self._system.set_property(traj_propname.source(), trajectory)

    def clear_energy_trajectory(self, map=None):
        """
        Completely clear any existing energy trajectory
        """
        from ..base import create_map

        map = create_map(map)

        traj_propname = map["energy_trajectory"]

        if traj_propname.has_value():
            return

        from ..maths import EnergyTrajectory

        self._system.set_property(traj_propname.source(), EnergyTrajectory())

    def evaluate(self, *args, **kwargs):
        """Return an evaluator for this Systme (or of the matching
        index/search subset of this System)"""
        return self.molecules().evaluate(*args, **kwargs)

    def has_property(self, *args, **kwargs):
        """Return whether or not this system has the passed property"""
        return self._system.contains_property(*args, **kwargs)

    def property(self, *args, **kwargs):
        """Return the System property that matches the passed key/arguments"""
        return self._system.property(*args, **kwargs)

    def set_property(self, *args, **kwargs):
        """Set the System property according to the passed arguments"""
        self._system.set_property(*args, **kwargs)
        self._molecules = None

    def set_shared_property(self, name, value):
        """Set the shared System property according to the passed arguments"""
        from ..base import wrap

        self._system.set_shared_property(name, wrap(value))
        self._molecules = None

    def add_shared_property(self, name, value=None):
        """
        Add the shared System property called 'name' to this system.
        If 'value' is supplied, then this also sets the property too.
        """
        if value is None:
            self._system.add_shared_property(name)
        else:
            from ..base import wrap

            self._system.add_shared_property(name, wrap(value))

        self._molecules = None

    def remove_shared_property(self, *args, **kwargs):
        """
        Completely remove the specified shared property, and set this
        as a non-shared property for the future
        """
        self._system.remove_shared_property(*args, **kwargs)
        self._molecules = None

    def remove_all_shared_properties(self, *args, **kwargs):
        """
        Completely remove all shared properties, and set no properties
        as being shared from this system
        """
        self._system.remove_all_shared_properties(*args, **kwargs)
        self._molecules = None

    def properties(self):
        """Return all of the System-level properties of this System"""
        return self._system.properties()

    def property_keys(self):
        """Return the keys (IDs) of all of the System-level properties
        of this System
        """
        return self._system.property_keys()

    def shared_properties(self):
        """
        Return all of the System properties that are being shared
        (copied) to all contained molecules
        """
        return self._system.shared_properties()

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

    def smarts(self, *args, **kwargs):
        """
        Return the molecule views in this container as smarts strings. Include
        hydrogens in 'include_hydrogens' is True. This returns a list
        of smarts strings, in the same order as the views in the container
        """
        return self.molecules().smarts(*args, **kwargs)

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
