from sire.legacy.Convert._SireOpenMM import EMLECallback


def test_callback():
    """Makes sure that a callback method works correctly"""

    class Test:
        def callback(self, a, b, c, d):
            return (sum(a + b + c + d), [a, b], [c, d])

    # Instantiate the class.
    test = Test()

    # Create a callback object.
    cb = EMLECallback(test, "callback")

    # Create some lists to hold test data.
    a = [1, 2, 3, 4]
    b = [5, 6, 7, 8]
    c = [9, 10, 11, 12]
    d = [13, 14, 15, 16]

    # Call the callback.
    result = cb.call(a, b, c, d)

    # Make sure the result is correct.
    assert (
        result
        == (136, [[1, 2, 3, 4], [5, 6, 7, 8]], [[9, 10, 11, 12], [13, 14, 15, 16]])
        == test.callback(a, b, c, d)
    )
