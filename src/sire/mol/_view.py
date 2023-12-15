__all__ = ["view"]

_nglview_import_error = None

try:
    import nglview as _nglview

    _has_nglview = True
except ImportError:
    _has_nglview = False
except AttributeError as e:
    _has_nglview = False
    _nglview_import_error = e


if _has_nglview:

    class _FrameCache:
        def __init__(self, max_frames=100):
            self._cache = {}
            self._order = []
            self._max_frames = min(0, max_frames)

        def __contains__(self, index):
            return index in self._cache

        def save(self, index, frame):
            while len(self._order) > self._max_frames:
                idx = self._order[0]
                del self._cache[idx]
                self._order = self._order[1:]

            if index not in self._cache:
                self._order.append(index)

            self._cache[index] = frame

        def is_empty(self):
            return len(self._order) == 0

        def clear(self):
            self._cache = {}
            self._order = {}

        def __getitem__(self, index):
            return self._cache[index]

    @_nglview.register_backend("sire")
    class _SireStructureTrajectory(_nglview.Trajectory, _nglview.Structure):
        def __init__(
            self,
            obj=None,
            align=None,
            frame=None,
            smooth=None,
            wrap=None,
            mapping=None,
            map=None,
        ):
            if type(obj) is _SireStructureTrajectory:
                self._traj = obj._traj
                self._map = obj._map
            elif obj is not None:
                from ._trajectory import TrajectoryIterator
                from ..system import System
                from ..base import create_map

                if System.is_system(obj):
                    # faster to view a SelectorMol
                    obj = obj.molecules()

                self._map = create_map(map)

                if type(obj) is TrajectoryIterator:
                    if align is True:
                        align = "*"
                    elif align is False:
                        align = None

                    if align is not None:
                        obj = obj.align(align=align, frame=frame)

                    if smooth not in [False, None, 0]:
                        obj = obj.smooth(smooth=smooth)

                    if wrap in [True, False]:
                        obj = obj.wrap(wrap)

                    self._traj = obj
                else:
                    self._traj = obj.trajectory(
                        align=align,
                        frame=frame,
                        smooth=smooth,
                        wrap=wrap,
                        mapping=mapping,
                        map=self._map,
                    )
            else:
                self._traj = None
                self._map = None

            # set 'coords_only' to True, so that the loading of
            # frames is quicker (prevents the slow assignment
            # of space and time to all molecules for each frame)
            self._map.set("coords_only", True)

            self.ext = "pdb"
            self.params = {}
            import uuid

            self.id = str(uuid.uuid4())
            self._cache = _FrameCache()

        def __repr__(self):
            return str(self)

        def __str__(self):
            # This is the additional string added by NGLView to
            # identify atoms and bonds.
            return ":"

        def get_structure_string(self):
            from .. import save_to_string

            return "\n".join(
                save_to_string(
                    self._traj.first(),
                    self.ext,
                    map={"use_atom_numbers": True},
                )
            )

        @property
        def n_frames(self):
            return max(1, len(self._traj))

        def get_coordinates(self, index):
            if index in self._cache:
                return self._cache[index]

            from ..io import get_coords_array
            from ..units import angstrom

            frame = self._traj[index].current()

            coords = get_coords_array(frame, units=angstrom, map=self._map)

            if coords.size == 0:
                return coords

            if self._cache.is_empty():
                # work out how many frames to cache based on size
                # (no more than 128 MB)
                max_frames = max(1, int(32 * 1024 * 1024 / coords.size))
                self._cache = _FrameCache(max_frames=max_frames)

            self._cache.save(index, coords)

            return coords

    class _Representations:
        def __init__(self, view):
            if type(view).__name__ == "TrajectoryIterator":
                view = view.current()

            self.view = view
            self.atoms = view.atoms()
            self.reps = {}
            self.supported_reps = [
                "ball_and_stick",
                "base",
                "cartoon",
                "distance",
                "hyperball",
                "licorice",
                "line",
                "none",
                "point",
                "ribbon",
                "rocket",
                "rope",
                "spacefill",
                "surface",
                "trace",
                "tube",
            ]

            self.rest = set(range(0, len(self.atoms)))

        def get_selection_string(self, selection):
            try:
                s = self.atoms.find(self.view[selection].atoms())
                s = [str(x) for x in s]
                return "@" + ",".join(s)
            except Exception as e:
                from ..utils import Console

                Console.warning(
                    f"Unrecognised selection '{selection}'\n" f"Error is {e}"
                )
                return "*"

        def add(self, selection, rep):
            if type(selection) is not list:
                selection = [selection]

            if rep is None:
                rep = "none"

            if type(rep) is not list:
                rep = [rep]

            for s in selection:
                if s == "all":
                    a = list(range(0, len(self.atoms)))
                else:
                    try:
                        atms = self.view[s]
                        a = self.atoms.find(atms.atoms())
                    except Exception:
                        continue

                for r in rep:
                    if r not in self.reps:
                        self.reps[r] = set()

                    self.reps[r] = self.reps[r].union(a)
                    self.rest = self.rest.difference(a)

        def _add_rep(self, view, typ, atoms):
            typ = typ.strip().lstrip().rstrip().lower()

            if len(typ) == 0:
                typ = "none"

            parts = typ.split(":")

            if len(parts) == 1:
                colour = None
                opacity = None
            elif len(parts) == 2:
                typ = parts[0]

                try:
                    opacity = float(parts[1])
                    colour = None
                except Exception:
                    opacity = None
                    colour = parts[1]
            else:
                typ = parts[0]

                try:
                    opacity = float(parts[1])
                    colour = parts[2]
                except Exception:
                    opacity = float(parts[2])
                    colour = parts[1]

            if typ not in self.supported_reps:
                s = ", ".join(self.supported_reps)

                raise KeyError(
                    f"Unsupported representation '{typ}'. "
                    f"Available representations are: [ {s} ]."
                )

            if typ == "none":
                # don't display anything
                return

            if typ == "ball_and_stick":
                typ = "ball+stick"

            atoms = list(atoms)

            if colour is None and opacity is None:
                view.add_representation(typ, selection=atoms)
            elif opacity is None:
                view.add_representation(typ, selection=atoms, color=colour)
            elif colour is None:
                view.add_representation(typ, selection=atoms, opacity=opacity)
            else:
                view.add_representation(
                    typ, selection=atoms, color=colour, opacity=opacity
                )

        def populate(self, view, rest=None):
            if rest is not None:
                if type(rest) is not list:
                    rest = [rest]

                for r in rest:
                    self._add_rep(view, r, self.rest)

            for key, value in self.reps.items():
                self._add_rep(view, key, value)

    def view(
        obj,
        default: bool = True,
        no_default: bool = False,
        orthographic: bool = True,
        protein: str = "",
        water: str = "",
        ions: str = "",
        rest: str = "",
        all: str = "",
        ball_and_stick: str = "",
        base: str = "",
        cartoon: str = "",
        hyperball: str = "",
        licorice: str = "",
        line: str = "",
        point: str = "",
        ribbon: str = "",
        rocket: str = "",
        rope: str = "",
        spacefill: str = "",
        surface: str = "",
        trace: str = "",
        tube: str = "",
        center: str = None,
        align: str = None,
        frame: int = None,
        smooth=False,
        wrap=True,
        mapping=None,
        stage_parameters: str = None,
        map=None,
    ):
        """
        Return an NGLView viewer for this view. The returned
        viewer can be passed directly to, e.g. a Jupyter notebook
        to directly view the molecule(s), or it can be captured
        in a variable so that it's NGLViewer member functions
        can be called to edit the viewer before display.

        Full instructions on how to use this view are in
        the cheatsheet (https://sire.openbiosim.org/cheatsheet/view)

        center:
          Pass in a selection string to select the atoms to center
          in the view. By default no atoms are centered

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
          Whether or not to wrap the coordinates into the periodic box.

        orthographic:
          Set to False to use a perspective view, or accept the default
          (True) for an orthographic view

        default:
          The default representation / color for the view. If this
          is set to False or None or "none" then default views
          are disabled.

        no_default:
          Set to False to disable all default views

        protein, water, ions:
          Set the default views for protein, water and ion molecules.
          Set to False, None or "none" to disable these default views.

        rest:
          Synonym for "default", but does not disable all default views
          if set to None, False or "none"

        all:
          Set the representation for all atoms. Set to None, False or "none"
          to disable all views.

        ball_and_stick, base, cartoon etc.
          Set the selection strings for atoms that should be represented
          with these views.

        stage_parameters: dict
          An optional dictionary that will be passed directly
          to the NGLView object to set the stage parameters.

        map: dict or sire.base.PropertyMap
          An optional property map that can be used to control
          which properties are used to get the molecular data
          to be viewed.
        """
        from ..utils import NullProfiler

        p = NullProfiler()
        p = p.start("total")

        p1 = p.start("Traj")
        struc_traj = _SireStructureTrajectory(
            obj,
            align=align,
            frame=frame,
            smooth=smooth,
            wrap=wrap,
            mapping=mapping,
            map=map,
        )
        p1.stop()

        p1 = p.start("widget")
        view = _nglview.NGLWidget(struc_traj)
        p1.stop()

        if orthographic:
            view.camera = "orthographic"
        else:
            view.camera = "perspective"

        view.clear_representations()

        p1 = p.start("representations")
        reps = _Representations(obj)
        p1.stop()

        p1 = p.start("check reps")

        if all is None or all is False or default is None or default is False:
            # User has turned off all representations
            no_default = True
            default = False

        elif len(all) > 0:
            # don't have any default representations if the user
            # has asked that all atoms have a particular
            # representation
            no_default = True
            default = False

        try:
            # Allow the user to use `default` to mean `rest`
            if len(default) != 0:
                rest = default
                default = True
        except Exception:
            pass

        if default and (not no_default):
            # Add the defaults here so that the user can
            # override them, and also so that we can
            # disable defaults if needed
            if (protein is not None) and (protein is not False):
                if len(protein) == 0:
                    protein = "cartoon:sstruc"

            if (water is not None) and (water is not False):
                if len(water) == 0:
                    water = "line:0.5"

            if (ions is not None) and (ions is not False):
                if len(ions) == 0:
                    # use this, as the spheres are smaller than spacefill
                    ions = "ball_and_stick"

            if (rest is not None) and (rest is not False):
                if len(rest) == 0:
                    rest = "hyperball"

        if (protein is not None) and (protein is not False):
            if len(protein) > 0:
                reps.add("protein", protein)
        else:
            reps.add("protein", None)

        if (water is not None) and (water is not False):
            if len(water) > 0:
                reps.add("water", water)
        else:
            reps.add("water", None)

        if (ions is not None) and (ions is not False):
            if len(ions) > 0:
                reps.add("molecules with count(atoms) == 1", ions)
        else:
            reps.add("ions", None)

        def _split(rep):
            parts = rep.split(":")

            if len(parts) == 1:
                return rep, ""
            else:
                return parts[0], ":" + ":".join(parts[1:])

        def _check(rep):
            if (rep is not None) and (rep is not False):
                if len(rep) > 0:
                    return rep

            return None

        ball_and_stick = _check(ball_and_stick)
        base = _check(base)
        cartoon = _check(cartoon)
        licorice = _check(licorice)
        hyperball = _check(hyperball)
        line = _check(line)
        point = _check(point)
        ribbon = _check(ribbon)
        rocket = _check(rocket)
        rope = _check(rope)
        spacefill = _check(spacefill)
        surface = _check(surface)
        trace = _check(trace)
        tube = _check(tube)
        all = _check(all)

        p1.stop()

        p2 = p.start("add representations")

        if ball_and_stick is not None:
            ball_and_stick, colours = _split(ball_and_stick)
            reps.add(ball_and_stick, f"ball_and_stick{colours}")

        if base is not None:
            base, colours = _split(base)
            reps.add(base, f"base{colours}")

        if cartoon is not None:
            cartoon, colours = _split(cartoon)
            reps.add(cartoon, f"cartoon{colours}")

        if licorice is not None:
            licorice, colours = _split(licorice)
            reps.add(licorice, f"licorice{colours}")

        if hyperball is not None:
            hyperball, colours = _split(hyperball)
            reps.add(hyperball, f"hyperball{colours}")

        if line is not None:
            line, colours = _split(line)
            reps.add(line, f"line{colours}")

        if point is not None:
            point, colours = _split(point)
            reps.add(point, f"point{colours}")

        if ribbon is not None:
            ribbon, colours = _split(ribbon)
            reps.add(ribbon, f"ribbon{colours}")

        if rocket is not None:
            rocket, colours = _split(rocket)
            reps.add(rocket, f"rocket{colours}")

        if rope is not None:
            rope, colours = _split(rope)
            reps.add(rope, f"rope{colours}")

        if spacefill is not None:
            spacefill, colours = _split(spacefill)
            reps.add(spacefill, f"spacefill{colours}")

        if surface is not None:
            surface, colours = _split(surface)
            reps.add(surface, f"surface{colours}")

        if trace is not None:
            trace, colours = _split(trace)
            reps.add(trace, f"trace{colours}")

        if tube is not None:
            tube, colours = _split(tube)
            reps.add(tube, f"tube{colours}")

        if all is not None:
            reps.add("all", all)

        p2.stop()

        p1 = p.start("stage parameters")

        if stage_parameters is None:
            view.stage.set_parameters(
                clipNear=0,
                clipFar=100,
                clipDist=0,
                fogNear=100,
                fogFar=1000,
                backgroundColor="black",
            )
        else:
            view.stage.set_parameters(**stage_parameters)

        if (rest is None) or (rest is False):
            rest = None

        p1.stop()

        p1 = p.start("populate reps")
        reps.populate(view, rest=rest)
        p1.stop()

        p1 = p.start("center")
        if center is not None:
            view.center(selection=reps.get_selection_string(center))
        else:
            view.center()
        p1.stop()

        p.stop()

        if not p.is_null():
            print(p)

        return view

elif _nglview_import_error is not None:

    def view(obj, *args, **kwargs):
        raise ImportError(
            "nglview cannot be imported. This is because of an error "
            f"when nglview was loaded ({_nglview_import_error}). One "
            "possibility is that nglview is incompatible with the installed "
            "version of ipywidgets. Try to downgrade ipywidgets, e.g. "
            "\"conda install 'ipywidgets>=7.6.0,<8'\". You will need to "
            "restart Python and run this script/notebook again."
        )

else:

    def view(obj, *args, **kwargs):
        raise ImportError(
            "You need to install nglview to be able to view "
            "molecules. Do this by typing, e.g. "
            "'conda install nglview' and then restarting Python "
            "and running this script/notebook again."
        )
