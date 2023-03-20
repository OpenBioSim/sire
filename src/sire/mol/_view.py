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
        def __init__(self, obj=None, map=None):
            if type(obj) is _SireStructureTrajectory:
                self._traj = obj._traj
                self._map = obj._map
            elif obj is not None:
                from ._trajectory import TrajectoryIterator
                from ..base import create_map

                self._map = create_map(map)

                if type(obj) is TrajectoryIterator:
                    self._traj = obj
                else:
                    self._traj = obj.trajectory(self._map)
            else:
                self._traj = None
                self._map = None

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

            return "\n".join(save_to_string(self._traj.first(), self.ext))

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
            self.view = view
            self.atoms = view.atoms()
            self.reps = {}
            self.supported_reps = [
                "ball_and_stick",
                "cartoon",
                "licorice",
                "line",
                "point",
                "ribbon",
                "rocket",
                "rope",
                "spacefill",
                "surface",
                "trace",
                "tube"
            ]

            self.rest = set(range(0, len(self.atoms)))

        def add(self, selection, rep):
            if type(selection) is not list:
                selection = [selection]

            for s in selection:
                try:
                    a = self.atoms.find(self.view[s].atoms())
                except Exception:
                    continue

                if rep not in self.reps:
                    self.reps[rep] = set()

                self.reps[rep] = self.reps[rep].union(a)
                self.rest = self.rest.difference(a)

        def _add_rep(self, view, typ, atoms):
            typ = typ.strip().lstrip().rstrip().lower()

            if typ not in self.supported_reps:
                s = ", ".join(self.supported_reps)

                raise KeyError(
                    f"Unsupported representation '{typ}'. "
                    f"Available representations are: [ {s} ]."
                )

            if typ == "ball_and_stick":
                typ = "ball+stick"

            atoms = list(atoms)

            view.add_representation(typ, selection=atoms)

        def populate(self, view, rest=None):
            if rest is not None:
                self._add_rep(view, rest, self.rest)

            for key, value in self.reps.items():
                self._add_rep(view, key, value)

    def view(obj,
             no_default: bool = False,
             orthographic: bool = True,
             protein:str ="cartoon",
             water:str="line",
             ions:str="spacefill",
             rest:str="licorice",
             all:str=None,
             ball_and_stick:str=None,
             cartoon:str=None,
             licorice:str=None,
             line:str=None,
             point:str=None,
             ribbon:str=None,
             rocket:str=None,
             rope:str=None,
             spacefill:str=None,
             surface:str=None,
             trace:str=None,
             tube:str=None,
             stage_parameters:str=None,
             map=None):
        """
        Return an NGLView viewer for this view. The returned
        viewer can be passed directly to, e.g. a Jupyter notebook
        to directly view the molecule(s), or it can be captured
        in a variable so that it's NGLViewer member functions
        can be called to edit the viewer before display.

        See the NGLView documentation for more information
        on how to configure the viewer.

        https://nglviewer.org/#nglview

         stage_parameters: dict
             An optional dictionary that will be passed directly
             to the NGLView object to set the stage parameters.

         map: dict or sire.base.PropertyMap
             An optional property map that can be used to control
             which properties are used to get the molecular data
             to be viewed.
        """
        struc_traj = _SireStructureTrajectory(obj, map=map)
        view = _nglview.NGLWidget(struc_traj)

        if orthographic:
            view.camera = "orthographic"
        else:
            view.camera = "perspective"

        view.clear_representations()

        atoms = obj.atoms()

        reps = _Representations(obj)

        if no_default:
            protein = None
            water = None
            ions = None
            rest = None

        if protein is not None:
            reps.add("protein", protein)

        if water is not None:
            reps.add("water", water)

        if ions is not None:
            reps.add("molecules with count(atoms) == 1", ions)

        if ball_and_stick is not None:
            reps.add(ball_and_stick, "ball_and_stick")

        if cartoon is not None:
            reps.add(cartoon, "cartoon")

        if licorice is not None:
            reps.add(licorice, "licorice")

        if line is not None:
            reps.add(line, "line")

        if point is not None:
            reps.add(point, "point")

        if ribbon is not None:
            reps.add(ribbon, "ribbon")

        if rope is not None:
            reps.add(rope, "rope")

        if spacefill is not None:
            reps.add(spacefill, "spacefill")

        if surface is not None:
            reps.add(surface, "surface")

        if trace is not None:
            reps.add(trace, "trace")

        if tube is not None:
            reps.add(tube, "tube")

        if all is not None:
            reps.add("all", all)

        if stage_parameters is None:
            view.stage.set_parameters(
                clipNear=0,
                clipFar=100,
                clipDist=0,
                fogNear=0,
                fogFar=100,
                backgroundColor="black",
            )
        else:
            view.stage.set_parameters(**stage_parameters)

        reps.populate(view, rest=rest)

        view.center()

        return view

elif _nglview_import_error is not None:

    def view(obj, *args, **kwargs):
        raise ImportError(
            "nglview cannot be imported. This is because of an error "
            f"when nglview was loaded ({_nglview_import_error}). One "
            "possibility is that nglview is incompatible with the installed "
            "version of ipywidgets. Try to downgrade ipywidgets, e.g. "
            "\"mamba install 'ipywidgets>=7.6.0,<8'\". You will need to "
            "restart Python and run this script/notebook again."
        )

else:

    def view(obj, *args, **kwargs):
        raise ImportError(
            "You need to install nglview to be able to view "
            "molecules. Do this by typing, e.g. "
            "'mamba install nglview' and then restarting Python "
            "and running this script/notebook again."
        )
