__all__ = ["PageCache"]

from ..legacy.Base import PageCache


def __cache__(obj, data):
    """
    Add the passed object onto the cache. This will convert the object
    into a binary form (pickled, then hex-encoded) and it will store
    it in the cache. This returns a handle to the object in the cache,
    which can be used to restore it.
    """
    from ..legacy.Qt import QByteArray
    from pickle import dumps

    data = dumps(data)

    # now convert this to a QByteArray
    b = QByteArray(data.hex())
    return obj.__orig__cache__(b)


def __fetch__(obj):
    """
    Fetch the object from the cache and return it
    """
    from ..legacy.Qt import QByteArray
    from pickle import loads

    data = obj.__orig__fetch__()

    if hasattr(data, "constData"):
        data = bytes.fromhex(data.constData())
    else:
        data = bytes.fromhex(data.const_data())

    return loads(data)


if not hasattr(PageCache, "__orig__cache__"):
    PageCache.__orig__cache__ = PageCache.cache
    PageCache.cache = __cache__

    if hasattr(PageCache, "handle"):
        PageCache.Handle = PageCache.handle
        delattr(PageCache, "handle")

    if hasattr(PageCache, "page"):
        PageCache.Page = PageCache.page
        delattr(PageCache, "page")

    PageCache.Handle.__orig__fetch__ = PageCache.Handle.fetch
    PageCache.Handle.fetch = __fetch__
