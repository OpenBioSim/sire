from functools import partialmethod

from sire.legacy.Convert._SireOpenMM import EMLECallback

# Create some lists to hold test data.
a = [1, 2, 3, 4]
b = [5, 6, 7, 8]
c = [9, 10, 11, 12]
d = [13, 14, 15, 16]


def callback(a, b, c, d):
    # Zip lists together and return the sum.
    result = []
    for aa, bb, cc, dd in zip(a, b, c, d):
        result.append(aa + bb + cc + dd)
    return result


def callback_wrapper(obj, a, b, c, d):
    """A callback wrapper function"""
    # No object, compute the result directly.
    if obj is None:
        return callback(a, b, c, d)
    # Object, call the method.
    else:
        return obj._callback(a, b, c, d)


def test_callback_function():
    """Makes sure that a callback function works correctly"""

    # Create a callback object.
    cb = EMLECallback(None, callback_wrapper)

    # Call the callback.
    result = cb.call(a, b, c, d)

    # Make sure the result is correct.
    assert (
        result
        == [28, 32, 36, 40]
        == callback_wrapper(None, a, b, c, d)
        == callback(a, b, c, d)
    )


def test_callback_method():
    """Makes sure that a callback method works correctly"""

    class Test:
        def _callback(self, a, b, c, d):
            return callback(a, b, c, d)

    # Instantiate the class.
    test = Test()

    # Create a callback object.
    cb = EMLECallback(test, callback_wrapper)

    # Call the callback.
    result = cb.call(a, b, c, d)

    # Make sure the result is correct.
    assert (
        result
        == [28, 32, 36, 40]
        == callback_wrapper(test, a, b, c, d)
        == callback(a, b, c, d)
    )
