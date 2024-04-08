__all__ = ["TrajectoryIterator"]

# make sure that GeneralUnit has been modernised
from ..units import GeneralUnit as _GeneralUnit  # noqa: F401


def _get_align_atoms_and_reference(
    view, align, frame=None, mapping=None, map=None, find_all=True
):
    if align is None:
        # this is an all-to-all alignment
        atoms = view.atoms()

        if frame is not None:
            atoms.load_frame(frame)

        return (atoms, atoms)

    elif hasattr(align, "atoms"):
        align = align.atoms()

        if mapping is None:
            from . import AtomMapping

            mapping = AtomMapping(align, align)

        # use the mapping to go from the atoms to align to
        # the atoms in the view
        atoms = mapping.find(align.atoms(), view.atoms(), find_all=find_all)

        if frame is not None:
            align.load_frame(frame)

        if len(atoms) != len(align) and not find_all:
            align = mapping.swap().find(atoms, align, find_all=True)

        return (atoms, align)

    else:
        # This is a search string in 'view'
        atoms = view[align].atoms()

        if frame is not None:
            atoms.load_frame(frame)

        return (atoms, atoms)


class TrajectoryIterator:
    """
    An iterator that can be used to control which frames of a trajectory
    are accessed or processed.
    """

    def __init__(
        self,
        view=None,
        align=None,
        smooth=None,
        wrap=None,
        mapping=None,
        frame=None,
        map=None,
    ):
        if view is not None:
            from ..base import create_map

            self._map = create_map(map)
            self._view = view
            self._values = list(range(0, max(1, self._view.num_frames(self._map))))
            self._times = None
            self._iter = None
            self._frame = None
            self._mapping = mapping
            self._frame = frame

            if align is True:
                # passing in 'True' means align against everything
                align = "*"
            elif align is False:
                # passing in 'False" means don't align anything
                align = None

            self._align = align

            # Get the values of smooth and wrap from either the
            # constructor args or the property map
            if (smooth is None) and self._map.specified("smooth"):
                smooth = self._map["smooth"].value().as_integer()

            if wrap is None and self._map.specified("wrap"):
                wrap = self._map["wrap"].value().as_boolean()

            if smooth not in [None, False, 0]:
                if smooth is True:
                    smooth = 5
                else:
                    smooth = int(smooth)

                    if smooth < 1:
                        smooth = 1
            else:
                smooth = 1

            if wrap in [None, False, 0]:
                wrap = False
            else:
                wrap = True

            # now place them back into the property map
            self._map.set("smooth", smooth)
            self._map.set("wrap", wrap)

            from ..legacy.Mol import TrajectoryAligner

            if align is not None:
                (align, reference) = _get_align_atoms_and_reference(
                    view=self._view,
                    align=align,
                    mapping=mapping,
                    frame=frame,
                    map=self._map,
                    find_all=False,
                )

                self._aligner = TrajectoryAligner(align, reference, map=self._map)
            elif wrap or (smooth != 1):
                self._aligner = TrajectoryAligner(
                    self._view.evaluate().center(),
                    map=self._map,
                )
            else:
                # We don't need to do anything to process the frames
                self._aligner = None
        else:
            self._view = None
            self._values = []
            self._times = None
            self._iter = None
            self._map = None
            self._frame = None
            self._align = None
            self._aligner = None
            self._mapping = None
            self._frame = None

    def __iter__(self):
        return self

    def __next__(self):
        if self._view is None or self._values is None:
            raise StopIteration()

        if self._iter is None:
            self._iter = self._values.__iter__()

        self._frame = self._iter.__next__()

        return self.current()

    def __len__(self):
        return len(self._values)

    def __getitem__(self, val):
        from copy import copy

        it = copy(self)

        if type(val) is int:
            it._values = [self._values[val]]
        elif type(val) is slice:
            it._values = self._values[val]
        else:
            it._values = [self._values[v] for v in val]

        return it

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if self._view is None:
            return "TrajectoryIterator::null"
        elif len(self._values) <= 1:
            return f"Trajectory({self._view}, num_frames=1)"
        else:
            return f"Trajectory({self._view}, num_frames={len(self._values)})"

    def _is_trajectory_iterator(self):
        """
        Return that this is a TrajectoryIterator
        """
        return True

    def _to_legacy_system(self):
        """
        Internal function used to convert the contained view into a
        legacy.Sire.System
        """
        from .._load import _to_legacy_system

        return _to_legacy_system(self.current())

    def frame_indexes(self):
        """
        Return the indexes of the frames in this trajectory that
        will be viewed
        """
        from copy import copy

        return copy(self._values)

    def align(self, align, frame=None):
        """
        Return a copy of this trajectory where each frame will be aligned
        against the atoms that match the search string 'align'
        """
        if frame is None:
            frame = self._frame

        t = TrajectoryIterator(
            view=self._view,
            align=align,
            frame=frame,
            mapping=self._mapping,
            map=self._map,
        )
        t._values = self._values
        return t

    def smooth(self, smooth):
        """
        Return a copy of this trajectory where each frame will be smoothed
        over the specified number of frames (or the recommended number
        if 'smooth' is set to 'True')
        """
        t = TrajectoryIterator(
            view=self._view,
            align=self._align,
            frame=self._frame,
            smooth=smooth,
            mapping=self._mapping,
            map=self._map,
        )
        t._values = self._values
        return t

    def wrap(self, autowrap=True):
        """
        Return a copy of this trajectory where each frame will be auto-wrapped
        into the current space
        """
        t = TrajectoryIterator(
            view=self._view,
            align=self._align,
            frame=self._frame,
            wrap=autowrap,
            mapping=self._mapping,
            map=self._map,
        )
        t._values = self._values
        return t

    def num_frames(self):
        return len(self._values)

    def current(self):
        """Return the current frame in the trajectory"""
        if self._view is None or self._values is None:
            raise StopIteration()

        if self._frame is None:
            self._frame = self._values[0]

        ret = self._view.clone()

        if self._aligner is None:
            map = self._map
        else:
            map = self._map.clone()
            map.set("transform", self._aligner[self._frame])

        ret.load_frame(self._frame, map=map)

        try:
            mol = ret.molecule()
        except Exception:
            mol = ret[0].molecule()

        time_property = self._map["time"]

        if mol.has_property(time_property):
            time = mol.property(time_property)
            ret.frame_time = lambda: time
        else:
            from ..units import picosecond

            ret.frame_time = lambda: 0 * picosecond

        ret.frame_index = lambda: self._frame

        return ret

    def _populate_map(self, map):
        """
        Internal function called by sire.save to populate the
        property map that will be used by MoleculeParser to
        extract each frame from this trajectory.

        Note that this populates a copy of the map, which is
        returned
        """
        from ..base import create_map

        # do it this way around, so that 'map' overwrites
        # anything in this trajectory's map
        map = create_map(self._map, map)

        # now add in the frames to save
        frames_to_write = self.frame_indexes()
        map.set("frames_to_write", frames_to_write)

        # and the aligner
        if self._aligner is not None:
            map.set("frame_aligner", self._aligner)

        return map

    def first(self):
        """Return the first frame in the trajectory"""
        if self._view is None or self._values is None:
            raise StopIteration()

        old_frame = self._frame

        self._frame = self._values[0]

        ret = self.current()

        self._frame = old_frame

        return ret

    def times(self):
        if self._times is not None:
            return self._times

        if self._view is None:
            return {}

        # load the times from the actual underlying trajectory data
        try:
            mol = self._view.molecule()
        except Exception:
            mol = self._view[0].molecule()

        traj = mol.property(self._map["trajectory"])

        self._times = []

        for idx in self._values:
            self._times.append(traj[idx].time())

        return self._times

    def energies(self, obj1=None, to_pandas=True, map=None):
        if self._view is None:
            return {}

        import numpy as np

        from ..system import create_forcefield
        from ..legacy.System import calculate_trajectory_energies
        from .._colname import colname
        from . import _to_molecules
        from ..base import create_map

        map = self._map.merge(create_map(map))

        colnames = []
        forcefields = []

        if obj1 is None:
            for v in self.first():
                colnames.append(colname(v))
                forcefields.append(create_forcefield(v, map=map))
        else:
            if type(obj1) is TrajectoryIterator:
                if obj1.num_frames() != self.num_frames():
                    raise ValueError(
                        "The two trajectories have a different "
                        "number of frames! "
                        f"{self.num_frames()} versus f{obj1.num_frames()}."
                    )

                obj1_mols = _to_molecules(obj1.first())

                for v in self.first():
                    colnames.append(colname(v))
                    forcefields.append(create_forcefield(v, obj1_mols, map=map))
            else:
                for v in self.first():
                    colnames.append(colname(v))
                    forcefields.append(
                        create_forcefield(v, _to_molecules(obj1), map=map)
                    )

        nframes = len(self)

        times = np.zeros(nframes, dtype=float)
        indexes = np.zeros(nframes, dtype=int)

        t = self.times()

        for i, idx in enumerate(self._values):
            times[i] = t[i].to_default()
            indexes[i] = idx

        time_unit = t[0].get_default().unit_string()
        energy_unit = None

        components = {}

        ff_nrgs = calculate_trajectory_energies(forcefields, self._values, map=map)

        for ff_idx in range(0, len(forcefields)):
            nrg = ff_nrgs[ff_idx][0]

            if ff_idx == 0:
                energy_unit = nrg.get_default().unit_string()

            components[colname(colnames[ff_idx], "total")] = np.zeros(
                nframes, dtype=float
            )

            for key in nrg.components().keys():
                components[colname(colnames[ff_idx], key)] = np.zeros(
                    nframes, dtype=float
                )

        for i in range(0, nframes):
            for ff_idx in range(0, len(forcefields)):
                nrg = ff_nrgs[ff_idx][i]
                components[colname(colnames[ff_idx], "total")][i] = nrg.to_default()

                for key in nrg.components().keys():
                    try:
                        components[colname(colnames[ff_idx], key)][i] = nrg[
                            key
                        ].to_default()
                    except KeyError:
                        k = colname(colnames[ff_idx], key)
                        components[k] = np.zeros(nframes, dtype=float)
                        components[k][i] = nrg[key].to_default()

        data = {}

        data["frame"] = indexes
        data["time"] = times

        colnames = list(components.keys())
        colnames.sort()

        for name in colnames:
            data[name] = components[name]

        if to_pandas:
            import pandas as pd

            df = pd.DataFrame(data)
            df.set_index("frame")

            df.time_unit = lambda: time_unit
            df.energy_unit = lambda: energy_unit

            def pretty_plot(x="time", y=None):
                if y is None:
                    y = colnames

                if x == "time":
                    xlabel = f"Time / {df.time_unit()}"
                elif x == "frame":
                    xlabel = "Frame"
                else:
                    xlabel = x

                ax = df.plot(x=x, y=y)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(f"Energy / {df.energy_unit()}")
                ax.legend(bbox_to_anchor=(1.05, 1.0))

            df.pretty_plot = pretty_plot

            return df

        return data

    def energy(self, obj1=None, to_pandas=True, map=None):
        if self._view is None:
            return {}

        import numpy as np

        from ..system import create_forcefield
        from ..legacy.System import calculate_trajectory_energy
        from ..base import create_map

        map = self._map.merge(create_map(map))

        if obj1 is None:
            ff = create_forcefield(self.first(), map=map)
        else:
            if type(obj1) is TrajectoryIterator:
                if obj1.num_frames() != self.num_frames():
                    raise ValueError(
                        "The two trajectories have a different "
                        "number of frames! "
                        f"{self.num_frames()} versus f{obj1.num_frames()}."
                    )

                ff = create_forcefield(self.first(), obj1.first(), map=map)
            else:
                ff = create_forcefield(self.first(), obj1, map=map)

        nframes = len(self)

        times = np.zeros(nframes, dtype=float)
        indexes = np.zeros(nframes, dtype=int)

        t = self.times()

        for i, idx in enumerate(self._values):
            times[i] = t[i].to_default()
            indexes[i] = idx

        time_unit = t[0].get_default().unit_string()
        energy_unit = None

        # calculate all the energies
        nrgs = calculate_trajectory_energy(ff, self._values, map)

        # convert the result into a pandas dataframe
        components = {}

        nrg = nrgs[0]
        energy_unit = nrg.get_default().unit_string()
        components["total"] = np.zeros(nframes, dtype=float)
        for key in nrg.components().keys():
            components[key] = np.zeros(nframes, dtype=float)

        for i in range(0, nframes):
            nrg = nrgs[i]
            components["total"][i] = nrg.to_default()

            for key in nrg.components().keys():
                try:
                    components[key][i] = nrg[key].to_default()
                except KeyError:
                    components[key] = np.zeros(nframes, dtype=float)
                    components[key][i] = nrg[key].to_default()

        data = {}

        data["frame"] = indexes
        data["time"] = times

        colnames = list(components.keys())
        colnames.sort()

        for colname in colnames:
            data[colname] = components[colname]

        if to_pandas:
            import pandas as pd

            df = pd.DataFrame(data)
            df.set_index("frame")

            df.time_unit = lambda: time_unit
            df.energy_unit = lambda: energy_unit

            def pretty_plot(x="time", y=None):
                if y is None:
                    y = colnames

                if x == "time":
                    xlabel = f"Time / {df.time_unit()}"
                elif x == "frame":
                    xlabel = "Frame"
                else:
                    xlabel = x

                ax = df.plot(x=x, y=y)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(f"Energy / {df.energy_unit()}")
                ax.legend(bbox_to_anchor=(1.05, 1.0))

            df.pretty_plot = pretty_plot

            return df

        return data

    def _simple_measures(self, to_pandas):
        from .._colname import colname

        if self._view is None:
            return {}

        uses_measures = None

        if hasattr(self._view, "measures"):
            uses_measures = True
        elif hasattr(self._view, "measure"):
            uses_measures = False
        else:
            raise AttributeError(
                f"This view ({self._view}) does not have a `.measure()` "
                "or `.measures()` function, so cannot be measured."
            )

        import numpy as np

        colnames = []
        columns = []

        nframes = len(self)

        times = np.zeros(nframes, dtype=float)
        indexes = np.zeros(nframes, dtype=int)

        time_unit = None
        measure_unit = None

        from ..base import ProgressBar

        if uses_measures:
            for view in self._view:
                colnames.append(colname(view))
                columns.append(np.zeros(nframes, dtype=float))

            with ProgressBar(total=nframes, text="Looping through frames") as progress:
                for idx, frame in enumerate(self.__iter__()):
                    for i, measure in enumerate(frame.measures(map=self._map)):
                        columns[i][idx] = measure.to_default()
                        times[idx] = frame.frame_time().to_default()

                        if measure_unit is None:
                            if not measure.is_zero():
                                measure_unit = measure.get_default().unit_string()

                        if time_unit is None:
                            time = frame.frame_time()
                            if not time.is_zero():
                                time_unit = time.get_default().unit_string()

                    indexes[idx] = frame.frame_index()
                    progress.set_progress(idx)
        else:
            colnames.append(colname(self._view))
            column = np.zeros(nframes, dtype=float)

            with ProgressBar(total=nframes, text="Looping through frames") as progress:
                for idx, frame in enumerate(self.__iter__()):
                    measure = frame.measure(map=self._map)
                    column[idx] = measure.to_default()
                    times[idx] = frame.frame_time().to_default()

                    if measure_unit is None:
                        if not measure.is_zero():
                            measure_unit = measure.get_default().unit_string()

                    if time_unit is None:
                        time = frame.frame_time()
                        if not time.is_zero():
                            time_unit = time.get_default().unit_string()

                    indexes[idx] = frame.frame_index()
                    progress.set_progress(idx)

            columns = [column]

        data = {}

        data["frame"] = indexes
        data["time"] = times

        for i in range(0, len(colnames)):
            data[colnames[i]] = columns[i]

        if to_pandas:
            import pandas as pd

            df = pd.DataFrame(data)
            df.set_index("frame")

            df.time_unit = lambda: time_unit
            df.measure_unit = lambda: measure_unit

            def pretty_plot(x="time", y=None):
                if y is None:
                    y = colnames

                if x == "time":
                    xlabel = f"Time / {df.time_unit()}"
                elif x == "frame":
                    xlabel = "Frame"
                else:
                    xlabel = x

                ax = df.plot(x=x, y=y)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(f"Size / {df.measure_unit()}")
                ax.legend(bbox_to_anchor=(1.05, 1.0))

            df.pretty_plot = pretty_plot

            return df

        return data

    def _custom_measures(self, func, to_pandas):
        if self._view is None:
            return {}

        if not type(func) is dict:
            func = {"custom": func}

        import numpy as np

        colnames = []
        columns = []

        nframes = len(self)

        times = np.zeros(nframes, dtype=float)
        indexes = np.zeros(nframes, dtype=int)

        time_unit = None
        measure_unit = None

        for key in func.keys():
            colnames.append(key)
            columns.append(np.zeros(nframes, dtype=float))

        from ..base import ProgressBar

        with ProgressBar(total=nframes, text="Looping through frames") as progress:
            for idx, frame in enumerate(self.__iter__()):
                for i, f in enumerate(func.values()):
                    measure = f(frame)
                    columns[i][idx] = measure.to_default()
                    times[idx] = frame.frame_time().to_default()

                    if measure_unit is None:
                        if not measure.is_zero():
                            measure_unit = measure.get_default().unit_string()

                    if time_unit is None:
                        time = frame.frame_time()
                        if not time.is_zero():
                            time_unit = time.get_default().unit_string()

                indexes[idx] = frame.frame_index()
                progress.set_progress(idx)

        data = {}

        data["frame"] = indexes
        data["time"] = times

        for i in range(0, len(colnames)):
            data[colnames[i]] = columns[i]

        if to_pandas:
            import pandas as pd

            df = pd.DataFrame(data)
            df.set_index("frame")

            df.time_unit = lambda: time_unit
            df.measure_unit = lambda: measure_unit

            def pretty_plot(x="time", y=None):
                if y is None:
                    y = colnames

                if x == "time":
                    xlabel = f"Time / {df.time_unit()}"
                elif x == "frame":
                    xlabel = "Frame"
                else:
                    xlabel = x

                ax = df.plot(x=x, y=y)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(f"Size / {df.measure_unit()}")
                ax.legend(bbox_to_anchor=(1.05, 1.0))

            df.pretty_plot = pretty_plot

            return df

        return data

    def measures(self, func=None, to_pandas=True):
        if func is None:
            return self._simple_measures(to_pandas=to_pandas)
        else:
            return self._custom_measures(func=func, to_pandas=to_pandas)

    def apply(self, func, *args, **kwargs):
        """
        Call the passed function on all frames of the trajectory,
        appending the result to a list of results, which
        is returned.

        The function can be either;

        1. a string containing the name of the function to call, or
        2. an actual function (either a normal function or a lambda expression)

        You can optionally pass in positional and keyword arguments
        here that will be passed to the function.

        Args:
            func (str or function): The function to be called, or the name
                                    of the function to be called.

        Returns:
            list: A list of the results, with one result per
                  frame in the trajectory
        """
        result = []

        from ..base import ProgressBar

        nframes = len(self)

        if str(func) == func:
            # we calling a named function
            with ProgressBar(total=nframes, text="Looping through frames") as progress:
                for i in range(0, nframes):
                    obj = self.__getitem__(i).current()
                    result.append(getattr(obj, func)(*args, **kwargs))
                    progress.set_progress(i + 1)

        else:
            # we have been passed the function to call
            with ProgressBar(total=nframes, text="Looping through frames") as progress:
                for i in range(0, nframes):
                    obj = self.__getitem__(i).current()
                    result.append(func(obj, *args, **kwargs))
                    progress.set_progress(i + 1)

        return result

    def view(self, *args, **kwargs):
        from ._view import view

        return view(self, *args, **kwargs)

    def rmsd(
        self,
        reference=None,
        align: bool = None,
        frame: int = None,
        mapping=None,
        match_all: bool = True,
        map=None,
        to_pandas: bool = False,
    ):
        """
        Calculate the RMSD across the frames of this trajectory, against
        either all atoms or, optionally, the passed view. You can specify
        whether to align the trajectory against the view
        (controlling which frame should be the source of alignment via `frame`)

        reference: molecule container or search string or index or range
           Anything that is valid to specify a sub-view of the trajectory,
           or that can be a molecular container itself. If this is a
           view, then the RMSD only to the atoms in the view is calculated.
           If the atoms aren't in this trajectory, then you need to provide
           a mapping that shows how you can map from the atoms in the
           view to the atoms in the trajectory.

        align: bool
           Whether or not to align the frames of the trajectory against
           the reference. This defaults to True

        frame: int
           The frame number of the reference against which to align the
           trajectory

        mapping: AtomMapping
           The mapping from atoms in the reference to atoms in the
           trajectory. This is only needed if the reference atoms
           are not contained in the atoms in this trajectory.

        match_all: bool
            Whether or not to find and match all of the atoms in
            'reference' to the atoms in this trajectory. Normally
            you do want to do this, so it defaults to True. If this
            is False, then the RMSD of only atoms in this trajectory
            that in 'reference' will be calculated, even if it just
            a single atom.

        map: dict
           Any parameters that will overwrite any of the map
           parameters that are already in this trajectory

        to_pandas: bool
           Whether or not to write the output to a Pandas DataFrame
           (default False)
        """
        from ..legacy.Mol import TrajectoryAligner
        from ..legacy.Mol import get_rmsd

        if reference is None:
            # nothing to align against
            align = False

        elif align is None:
            # auto-align only for smaller systems
            if hasattr(reference, "molecules") and len(reference.molecules()) > 100:
                align = False
            else:
                align = True

        (atoms, reference) = _get_align_atoms_and_reference(
            view=self.current(),
            align=reference,
            frame=frame,
            mapping=mapping,
            map=self._map,
            find_all=match_all,
        )

        if map is None:
            map = self._map
        else:
            from ..base import create_map

            map = self._map.merge(create_map(map))

        if align:
            aligner = TrajectoryAligner(atoms, reference, map=map)
        else:
            aligner = TrajectoryAligner()

        rmsd = get_rmsd(
            atoms=atoms,
            reference=reference,
            aligner=aligner,
            frames=self._values,
            map=map,
        )

        if to_pandas:
            import pandas as pd

            if len(rmsd) == 0:
                return pd.DataFrame()

            # convert the RMSD to default length units
            from ..units import angstrom, picosecond

            length_unit = angstrom.get_default()
            time_unit = picosecond.get_default()

            rmsd = [x.to(length_unit) for x in rmsd]

            # get the times for each frame
            import copy

            t = copy.copy(self)
            t._view = self._view.atoms()[0]

            times = []

            for frame in t:
                times.append(frame.frame_time())

            times = [x.to(time_unit) for x in times]

            df = pd.DataFrame(
                {"frame": copy.copy(self._values), "time": times, "rmsd": rmsd}
            )

            df.time_unit = lambda: time_unit.unit_string()
            df.measure_unit = lambda: length_unit.unit_string()

            def pretty_plot(x="time", y=None):
                if y is None:
                    y = ["rmsd"]

                if x == "time":
                    xlabel = f"Time / {df.time_unit()}"
                elif x == "frame":
                    xlabel = "Frame"
                else:
                    xlabel = x

                ax = df.plot(x=x, y=y)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(f"RMSD / {df.measure_unit()}")
                ax.get_legend().remove()

            df.pretty_plot = pretty_plot

            return df
        else:
            return rmsd
