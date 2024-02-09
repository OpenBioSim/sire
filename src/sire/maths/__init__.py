__all__ = [
    "align",
    "create_quaternion",
    "get_alignment",
    "EnergyTrajectory",
    "kabasch",
    "kabasch_fit",
    "Matrix",
    "pi",
    "rotate",
    "RanGenerator",
    "Sphere",
    "Torsion",
    "Transform",
    "Triangle",
    "Vector",
]

from ..legacy import Maths as _Maths
from ..legacy.Maths import (
    Matrix,
    Quaternion,
    RanGenerator,
    Triangle,
    Transform,
    Torsion,
    pi,
    EnergyTrajectory,
)

from ._vector import Vector
from ._sphere import Sphere

from .. import use_new_api as _use_new_api

_use_new_api()


try:
    kabasch_fit = _Maths.kabaschFit
    get_alignment = _Maths.getAlignment
except AttributeError:
    kabasch_fit = _Maths.kabasch_fit
    get_alignment = _Maths.get_alignment


kabasch = _Maths.kabasch
align = _Maths.align


def rotate(point, angle=None, axis=None, matrix=None, quaternion=None, center=None):
    """
    Rotate the passed point by the passed angle and axis, or the
    passed matrix, or the passed quaternion, optionally centering
    the rotation about the passed center.

     point: sire.maths.Vector
         The vector to be rotated

     angle: (float or angle)
         The angle to rotate by - this is interpreted as
         degrees if you pass in a float. Otherwise use
         sire.units.degrees or sire.units.radians to specify
         the angle unit. This is superseded by the
         matrix and quaternion arguments.

     axis: sire.maths.Vector (or anything that can convert to a Vector)
         The vector about which to rotate. If this is not
         specified, and no other rotation specification is
         used, then the rotation is about the z axis.
         This is superseded by the matrix and
         quaternion arguments.

     quaternion: sire.maths.Quaternion
         The Quaternion description of the rotation. Note that,
         if you pass this, then the angle, axis and matrix
         arguments will be ignored.

     matrix: sire.maths.Matrix
         The 3x3 rotation matrix that describes the rotation.
         Note that, if you pass this, then the angle and axis
         arguments will be ignored. This is superseded by
         the quaternion argument.

     center: sire.maths.Vector
         The center of rotation. Defaults to (0,0,0) if not set.

     Returns: sire.maths.Vector
         The rotated vector
    """
    q = create_quaternion(angle=angle, axis=axis, matrix=matrix, quaternion=quaternion)

    if center is None:
        center = Vector(0, 0, 0)

    return _Maths.rotate(point, q.to_matrix(), center)


def create_quaternion(angle=None, axis=None, matrix=None, quaternion=None):
    """Create a quaternion from the passed angle and axis
    of the passed rotation matrix. If a rotation
    matrix is passed then this will ignore the
    passed angle and axis. If a quaternion is passed
    then this will ignore the matrix, angle and axis
    arguments.

     angle: (float or angle)
         The angle to rotate by - this is interpreted as
         degrees if you pass in a float. Otherwise use
         sire.units.degrees or sire.units.radians to specify
         the angle unit. This is superseded by the
         matrix and quaternion arguments.

     axis: sire.maths.Vector (or anything that can convert to a Vector)
         The vector about which to rotate. If this is not
         specified, and no other rotation specification is
         used, then the rotation is about the z axis.
         This is superseded by the matrix and
         quaternion arguments.

     quaternion: sire.maths.Quaternion
         The Quaternion description of the rotation. Note that,
         if you pass this, then the angle, axis and matrix
         arguments will be ignored.

     matrix: sire.maths.Matrix
         The 3x3 rotation matrix that describes the rotation.
         Note that, if you pass this, then the angle and axis
         arguments will be ignored. This is superseded by
         the quaternion argument.

     Returns: sire.maths.Quaternion
         The quaternion that represents the rotation
    """
    if quaternion is None:
        if type(angle) is Quaternion:
            # the user has passed in a quaternion as the first argument
            return angle

        if type(angle) is Matrix and matrix is None:
            # the user has passed in a rotation matrix as the first argument
            matrix = angle
            angle = None
            axis = None

        if matrix is None:
            if angle is None:
                raise ValueError(
                    "You must specify either the angle, rotation matrix "
                    "or quaternion used to rotate the molecule."
                )

            from ..units import degrees

            try:
                valid_angle = angle.has_same_units(degrees)
            except Exception:
                try:
                    angle = float(angle) * degrees
                    valid_angle = True
                except Exception:
                    valid_angle = False

            if not valid_angle:
                raise TypeError(
                    f"The passed angle of rotation ({angle}) has the wrong "
                    f"type ({type(angle)}). It should be an angle or a float."
                )

            if axis is None:
                axis = Vector(0, 0, 1)

            # construct from the passed angle and vector
            return Quaternion(angle, Vector(axis))
        else:
            if angle is not None or axis is not None:
                from ..utils import Console

                Console.warning(
                    "The angle and/or axis of rotation will be ignored "
                    "because you have passed in a rotation matrix."
                )

            if type(matrix) is not Matrix:
                raise TypeError(
                    f"The rotation matrix ({matrix}) must be of type "
                    f"sire.maths.Matrix. Type {type(matrix)} is not "
                    "supported."
                )

            return Quaternion(matrix)
    else:
        if matrix is not None:
            from ..utils import Console

            Console.warning(
                "The rotation matrix will be ignored "
                "because you have passed in a quaternion."
            )

        if angle is not None or axis is not None:
            from ..utils import Console

            Console.warning(
                "The angle and/or axis of rotation will be ignored "
                "because you have passed in a quaternion."
            )

        if type(quaternion) is not Quaternion:
            raise TypeError(
                f"The quaternion ({quaternion}) must be of type "
                f"sire.maths.Quaternion. Type {type(quaternion)} is not "
                "supported."
            )

        return quaternion


if not hasattr(EnergyTrajectory, "to_pandas"):

    def _to_pandas(
        obj, temperature=None, to_alchemlyb: bool = False, energy_unit: str = None
    ):
        """
        Return the energy trajectory as a pandas DataFrame

        Parameters
        ----------

        temperature: temperature
            The temperature of the simulation. If this is
            not set then the temperature from this table's
            `ensemble` or `temperature` property will be
            used. Note that you only need a temperature
            if you are converting to alchemlyb format.

        to_alchemlyb: bool
            This will format the DataFrame in a way that is
            compatible with alchemlyb. This will allow the
            DataFrame to be used as part of an alchemlyb
            free energy calculation.

        energy_unit: str
            Whichever of the alchemlyb energy units you want the output
            DataFrame to use. This is in alchemlyb format, e.g.
            `kcal/mol`, `kJ/mol`, or `kT`. This is only used if
            `to_alchemlyb` is set to True.
        """
        import pandas as pd
        from ..units import picosecond, kcal_per_mol, kJ_per_mol

        data = {}

        convert_to_kt = False

        if to_alchemlyb:
            time_unit = picosecond
            time_unit_string = "ps"

            if energy_unit == "kT":
                energy_unit = kcal_per_mol
                energy_unit_string = "kT"
                convert_to_kt = True
            elif energy_unit == "kJ/mol":
                energy_unit = kJ_per_mol
                energy_unit_string = "kJ/mol"
            else:
                energy_unit = kcal_per_mol
                energy_unit_string = "kcal/mol"

            if temperature is None:
                # look for the temperature in the ensemble property
                if obj.has_property("ensemble"):
                    temperature = obj.property("ensemble").temperature()

                # ok, try the temperature property
                if temperature is None and obj.has_property("temperature"):
                    temperature = obj.property("temperature")

                if temperature is None:
                    raise ValueError(
                        "You must specify the temperature of the simulation "
                        "when converting to alchemlyb format, or ensure that "
                        "the trajectory has an ensemble or temperature "
                        "property."
                    )
        else:
            time_unit = picosecond.get_default()
            time_unit_string = time_unit.unit_string()

            energy_unit = kcal_per_mol.get_default()
            energy_unit_string = energy_unit.unit_string()

        data["time"] = obj.times(time_unit)

        keys = obj.label_keys()
        keys.sort()

        if convert_to_kt:
            from ..units import k_boltz

            kT = (k_boltz * temperature).to(kcal_per_mol)

        for key in keys:
            if to_alchemlyb and key == "lambda":
                data["fep-lambda"] = obj.labels_as_numbers(key)
            else:
                # use float keys if possible
                try:
                    column_header = float(key)
                except Exception:
                    column_header = key

                try:
                    data[column_header] = obj.labels_as_numbers(key)
                except Exception:
                    data[column_header] = obj.labels(key)

        keys = obj.keys()
        keys.sort()

        if to_alchemlyb:
            keys.remove("kinetic")
            keys.remove("potential")

        for key in keys:
            # use float keys if possible
            try:
                column_header = float(key)
            except Exception:
                column_header = key

            nrgs = obj.energies(key, energy_unit)

            if convert_to_kt:
                nrgs = [x / kT for x in nrgs]

            data[column_header] = nrgs

        if to_alchemlyb:
            df = pd.DataFrame(data).set_index(["time", "fep-lambda"])
        else:
            df = pd.DataFrame(data).set_index("time")

        if temperature is not None:
            from .. import u
            from ..units import kelvin

            df.attrs["temperature"] = u(temperature).to(kelvin)

        df.attrs["energy_unit"] = energy_unit_string
        df.attrs["time_unit"] = time_unit_string

        return df

    def _to_alchemlyb(obj, temperature=None, energy_unit: str = "kcal/mol"):
        """
        Return the energy trajectory as an alchemlyb-formatted pandas DataFrame

        Parameters
        ----------

        temperature: temperature
            The temperature of the simulation. If this is
            not set then the temperature from this table's
            `ensemble` or `temperature` property will be
            used.

        energy_unit: str
            Whichever of the alchemlyb energy units you want the output
            DataFrame to use. This is in alchemlyb format, e.g.
            `kcal/mol`, `kJ/mol`, or `kT`

        Returns
        -------

        pandas.DataFrame
            A pandas DataFrame that is compatible with alchemlyb.
        """
        return obj.to_pandas(
            temperature=temperature, to_alchemlyb=True, energy_unit=energy_unit
        )

    EnergyTrajectory.to_pandas = _to_pandas
    EnergyTrajectory.to_alchemlyb = _to_alchemlyb
