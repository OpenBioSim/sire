from sire.legacy.Convert._SireOpenMM import EMLECallback


def test_callback():
    """Makes sure that a callback method works correctly"""

    class Test:
        def callback(self, a, b, c, d):
            return (42, d, c)

    # Instantiate the class.
    test = Test()

    # Create a callback object.
    cb = EMLECallback(test, "callback")

    # Create some lists to hold test data.
    a = [1, 2]
    b = [3, 4]
    c = [a, b]
    d = [b, a]

    # Call the callback.
    result = cb.call(a, b, c, d)

    # Make sure the result is correct.
    assert result == (42, d, c) == test.callback(a, b, c, d)
