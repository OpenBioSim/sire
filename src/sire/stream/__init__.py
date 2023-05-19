__all__ = [
    "save",
    "load",
    "get_data_header",
    "set_header_property",
    "get_header_property",
]

from ..legacy import Stream as _Stream

from .. import use_new_api as _use_new_api

_use_new_api()


def save(obj, filename=None):
    """
    Save the passed object to the sire streamed data format.
    If 'filename' is passed then the data is written to this file.
    Otherwise the data is returned as a binary array.
    """
    if filename is None:
        return _Stream.save(obj)
    else:
        _Stream.save(obj, filename)


def load(data):
    """
    Load the passed data from the sire streamed data format.
    If 'data' is a string, then this will load the appropriate
    file. Otherwise, it will assume data is a binary array,
    so will load the data directly from that.
    """
    return _Stream.load(data)


def get_data_header(data):
    """
    Return the header data from the sire streamed data.
    If 'data' is a string then this is loaded from the
    appropriate file. Otherwise it will assume data is a
    binary array and will load directly from that.
    """
    try:
        return _Stream.getDataHeader(data)
    except AttributeError:
        pass

    return _Stream.get_data_header(data)


def _to_binary(value):
    """Internal function to creata a binary array from 'value'"""
    import pickle

    return pickle.dumps(value).hex()


def _from_binary(data):
    """Internal function to create a value from the passed binary data"""
    import pickle

    return pickle.loads(bytes.fromhex(data))


def set_header_property(key, value):
    """
    Set the global header value for key 'key' to the passed 'value'.
    This will write this data into all sire streamed data files that
    are written from now on.
    """
    try:
        _Stream.setHeaderProperty(key, _to_binary(value))
        return
    except AttributeError:
        pass

    _Stream.set_header_property(key, _to_binary(value))


def get_header_property(key):
    """
    Return the global header value for the key 'key'. Returns
    None if no value is set.
    """
    try:
        return _from_binary(_Stream.getHeaderProperty(key))
    except AttributeError:
        pass

    return _from_binary(_Stream.get_header_property(key))


FileHeader = _Stream.FileHeader


def _fix_file_header():
    FileHeader.__orig__property = FileHeader.property

    def __get_property__(obj, key, default_value=None):
        try:
            data = obj.__orig__property(key)
        except Exception:
            return default_value

        return _from_binary(data)

    FileHeader.property = __get_property__


if not hasattr(FileHeader, "__orig__property"):
    _fix_file_header()
