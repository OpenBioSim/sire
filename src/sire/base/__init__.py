__all__ = [
    "create_map",
    "wrap",
    "PageCache",
    "PropertyMap",
    "ProgressBar",
    "Properties",
]

from ..legacy import Base as _Base

from .. import use_new_api as _use_new_api

from ..legacy.Base import PropertyMap, Properties

from ._progressbar import ProgressBar

from ._pagecache import PageCache

_use_new_api(is_base=True)

wrap = _Base.wrap


def _fix_wrapped_property_base_type(CLASS):
    """Private function that adds in all of the operators so that,
    as much as possible, the base type PropertyWrappers (e.g.
    NumberProperty, StringProperty) behave like numbers or
    strings, and can auto-convert to those types as needed.
    """

    def __value(v):
        try:
            v = v.value()
        except Exception:
            pass

        try:
            if v.is_integer():
                return int(v)
        except Exception:
            pass

        return v

    def __add__(obj, other):
        return __value(obj) + __value(other)

    CLASS.__add__ = __add__
    CLASS.__radd__ = __add__

    def __sub__(obj, other):
        return __value(obj) - __value(other)

    def __rsub__(obj, other):
        return __value(other) - __value(obj)

    CLASS.__sub__ = __sub__
    CLASS.__rsub__ = __rsub__

    def __eq__(obj, other):
        return __value(obj) == __value(other)

    CLASS.__eq__ = __eq__

    def __ne__(obj, other):
        return __value(obj) != __value(other)

    CLASS.__ne__ = __ne__

    def __gt__(obj, other):
        return __value(obj) > __value(other)

    CLASS.__gt__ = __gt__

    def __ge__(obj, other):
        return __value(obj) >= __value(other)

    CLASS.__ge__ = __ge__

    def __lt__(obj, other):
        return __value(obj) < __value(other)

    CLASS.__lt__ = __lt__

    def __le__(obj, other):
        return __value(obj) <= __value(other)

    CLASS.__le__ = __le__

    def __float__(obj):
        return float(obj.value())

    CLASS.__float__ = __float__

    def __int__(obj):
        return int(obj.value())

    CLASS.__int__ = __int__

    def __str__(obj):
        return str(__value(obj))

    CLASS.__str__ = __str__


_fix_wrapped_property_base_type(_Base.BooleanProperty)
_fix_wrapped_property_base_type(_Base.NumberProperty)
_fix_wrapped_property_base_type(_Base.StringProperty)


if not hasattr(PropertyMap, "__orig__set"):
    PropertyMap.__str__ = lambda x: str(x.to_dict())
    PropertyMap.__repr__ = PropertyMap.__str__

    def __propertymap_set(obj, key, value):
        try:
            obj.__orig__set(key, value)
            return
        except Exception:
            pass

        try:
            obj.__orig__set(key, _Base.PropertyName(value))
            return
        except Exception:
            pass

        obj.__orig__set(key, _Base.PropertyName(wrap(value)))

    PropertyMap.__orig__set = PropertyMap.set
    PropertyMap.set = __propertymap_set

    def __propertymap_get_string(obj, key: str):
        """
        Return the string value associated with the passed 'key'

        This returns 'key' if there is no value associated
        """
        key = str(key)
        if obj.specified(key):
            val = obj[key]

            if val.has_value():
                return val.value().as_string()
            else:
                return val.source()
        else:
            return key

    PropertyMap.get_string = __propertymap_get_string


def create_map(values, extras=None):
    """Construct a PropertyMap from the
    passed values. A PropertyMap is a class that lets you either provide
    extra options to some of the C++ functions, or to
    map the default keys used to find properties to
    your own keys

    You normally wouldn't use the class yourself.
    Instead, objects of this class will be created automatically
    from dictionaries, e.g.

    >>> mol.energy(map={"cutoff": 5*angstrom})

    would automatically create a PropertyMap, and is
    equivalent to writing

    >>> mol.energy(map=create_map({"cutoff": 5*angstrom}))

    In the above case you are providing an extra "cutoff"
    option, and are passing the value "5*angstrom".

    You can also use the map to change the keys used to
    find properties, e.g.

    >>> mol.energy(map={"coordinates": "my_coords"})

    would map the default "coordinates" property to your
    "my_coords" property. This means that the coordinates
    for the energy would be found at "my_coords" rather
    than "coordinates".

    You can map as many properties, and provide as many
    extra options as you want.
    """
    if values is None and extras is None:
        return PropertyMap()
    elif values is None:
        return create_map(extras)

    try:
        map = PropertyMap(values)
    except Exception:
        map = None

    if map is None:
        wrapped_values = {}

        for key, value in values.items():
            try:
                wrapped_values[key] = _Base.PropertyName(value)
            except Exception:
                try:
                    wrapped_values[key] = _Base.PropertyName(wrap(value))
                except Exception:
                    raise ValueError(
                        f"Unable to set the property {key} to {value} "
                        "because we cannot convert this to a C++ type."
                    )

        map = PropertyMap(wrapped_values)

    if extras is not None:
        return map + create_map(extras)
    else:
        return map
