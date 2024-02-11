from ..legacy.Maths import Vector as _Vector


def _fix_vector():
    def x(obj):
        from ..units import angstrom

        return obj.__old__x() * angstrom

    _Vector.__old__x = _Vector.x
    _Vector.x = x

    def y(obj):
        from ..units import angstrom

        return obj.__old__y() * angstrom

    _Vector.__old__y = _Vector.y
    _Vector.y = y

    def z(obj):
        from ..units import angstrom

        return obj.__old__z() * angstrom

    _Vector.__old__z = _Vector.z
    _Vector.z = z

    def __getitem__(obj, i):
        from ..units import angstrom

        return obj.__old__getitem__(i) * angstrom

    _Vector.__old__getitem__ = _Vector.__getitem__
    _Vector.__getitem__ = __getitem__

    def at(obj, i):
        from ..units import angstrom

        return obj.__old__at(i) * angstrom

    _Vector.__old_at = _Vector.at
    _Vector.at = at

    def getitem(obj, i):
        from ..units import angstrom

        return obj.__old__getitem(i) * angstrom

    _Vector.__old_getitem = _Vector.getitem
    _Vector.getitem = getitem

    def manhattan_length(obj):
        from ..units import angstrom

        return obj.__old_manhattan_length() * angstrom

    try:
        _Vector.__old_manhattan_length = _Vector.manhattan_length
    except Exception:
        _Vector.__old_manhattan_length = _Vector.manhattanLength
        delattr(_Vector, "manhattanLength")

    _Vector.manhattan_length = manhattan_length

    def length(obj):
        from ..units import angstrom

        return obj.__old_length() * angstrom

    _Vector.__old_length = _Vector.length
    _Vector.length = length

    def length2(obj):
        from ..units import angstrom

        return obj.__old_length2() * angstrom * angstrom

    _Vector.__old_length2 = _Vector.length2
    _Vector.length2 = length2

    def inv_length(obj):
        from ..units import angstrom

        return obj.__old_inv_length() / angstrom

    try:
        _Vector.__old_inv_length = _Vector.inv_length
    except Exception:
        _Vector.__old_inv_length = _Vector.invLength
        delattr(_Vector, "invLength")

    _Vector.inv_length = inv_length

    def inv_length2(obj):
        from ..units import angstrom

        return obj.__old_inv_length2() / (angstrom * angstrom)

    try:
        _Vector.__old_inv_length2 = _Vector.inv_length2
    except Exception:
        _Vector.__old_inv_length2 = _Vector.invLength2
        delattr(_Vector, "invLength2")

    _Vector.inv_length2 = inv_length2

    def distance(v1, v2):
        from ..units import angstrom

        v1 = Vector.to_vector(v1)
        v2 = Vector.to_vector(v2)

        return _Vector.__old_distance(v1, v2) * angstrom

    _Vector.__old_distance = _Vector.distance
    _Vector.distance = distance

    def distance2(v1, v2):
        from ..units import angstrom

        v1 = Vector.to_vector(v1)
        v2 = Vector.to_vector(v2)

        return _Vector.__old_distance2(v1, v2) * angstrom * angstrom

    _Vector.__old_distance2 = _Vector.distance2
    _Vector.distance2 = distance2

    def inv_distance(v1, v2):
        from ..units import angstrom

        v1 = Vector.to_vector(v1)
        v2 = Vector.to_vector(v2)

        return _Vector.__old_inv_distance(v1, v2) / angstrom

    try:
        _Vector.__old_inv_distance = _Vector.inv_distance
    except Exception:
        _Vector.__old_inv_distance = _Vector.invDistance
        delattr(_Vector, "invDistance")

    _Vector.inv_distance = inv_distance

    def inv_distance2(v1, v2):
        from ..units import angstrom

        v1 = Vector.to_vector(v1)
        v2 = Vector.to_vector(v2)

        return _Vector.__old_inv_distance2(v1, v2) / (angstrom * angstrom)

    try:
        _Vector.__old_inv_distance2 = _Vector.inv_distance2
    except Exception:
        _Vector.__old_inv_distance2 = _Vector.invDistance2
        delattr(_Vector, "invDistance2")

    _Vector.inv_distance2 = inv_distance2

    def __str__(obj):
        return f"( {obj.x()}, {obj.y()}, {obj.z()} )"

    _Vector.__str__ = __str__
    _Vector.__repr__ = __str__

    def __format__(obj, format):
        return (
            f"({obj.x().value():{format}}, {obj.y().value():{format}}, "
            f"{obj.z().value():{format}})"
        )

    _Vector.__format__ = __format__

    def __add__(obj, other):
        try:
            other = Vector(other)
        except Exception:
            pass

        return obj.__old__add__(other)

    _Vector.__old__add__ = _Vector.__add__
    _Vector.__add__ = __add__

    def __sub__(obj, other):
        try:
            other = Vector(other)
        except Exception:
            pass

        return obj.__old__sub__(other)

    def __rsub__(obj, other):
        try:
            other = Vector(other)
        except Exception:
            pass

        return other.__old__sub__(obj)

    _Vector.__old__sub__ = _Vector.__sub__
    _Vector.__sub__ = __sub__

    _Vector.__radd__ = __add__
    _Vector.__rsub__ = __rsub__

    def __iadd__(obj, other):
        try:
            other = Vector(other)
        except Exception:
            pass

        return obj.__old__iadd__(other)

    _Vector.__old__iadd__ = _Vector.__iadd__
    _Vector.__iadd__ = __iadd__

    def __isub__(obj, other):
        try:
            other = Vector(other)
        except Exception:
            pass

        return obj.__old__isub__(other)

    _Vector.__old__isub__ = _Vector.__isub__
    _Vector.__isub__ = __isub__


class Vector(_Vector):
    """
    A 3D point in space, or a 3D vector in space. This is a simple
    class containing 3 double precision values. These values
    represent points in units of Angstroms.
    """

    def __init__(self, *args, **kwargs):
        from ..units import angstrom
        from .. import u

        # get the default unit of length
        length_unit = angstrom.get_default()

        def _is_number(v):
            return isinstance(v, int) or isinstance(v, float)

        def _is_str(v):
            return isinstance(v, str)

        # mix of doubles and lengths?
        new_args = []

        for i in range(0, len(args)):
            if _is_number(args[i]):
                new_args.append((args[i] * length_unit).to(angstrom))
            elif _is_str(args[i]):
                new_args.append(u(args[i]).to(angstrom))
            elif hasattr(args[i], "__len__"):
                new_arg = []
                for arg in args[i]:
                    if _is_number(arg):
                        new_arg.append((arg * length_unit).to(angstrom))
                    elif _is_str(arg):
                        new_arg.append(u(arg).to(angstrom))
                    else:
                        new_arg.append(arg.to(angstrom))
                new_args.append(new_arg)
            else:
                new_args.append(args[i].to(angstrom))

        for key in kwargs.keys():
            if _is_number(kwargs[key]):
                kwargs[key] = (kwargs[key] * length_unit).to(angstrom)
            elif _is_str(kwargs[key]):
                kwargs[key] = u(kwargs[key]).to(angstrom)
            else:
                kwargs[key] = kwargs[key].to(angstrom)

        super().__init__(*new_args, **kwargs)

    @staticmethod
    def to_vector(arg):
        if isinstance(arg, _Vector) or isinstance(arg, Vector):
            return arg

        try:
            return Vector(arg)
        except Exception:
            pass

        try:
            return Vector(*arg)
        except Exception:
            pass

        raise TypeError(f"Could not convert {arg} to a sire.maths.Vector!")


if not hasattr(_Vector, "__old__getitem__"):
    _fix_vector()
