__all__ = ["Option"]

from enum import Enum as _Enum


class Option(str, _Enum):
    """
    Base class of all of the Option objects. These are
    effectively StrEnum objects that have additional
    functions to support case-insentive assigment
    and querying.

    It is very easy to create a new type of Option.

    Example
    -------

    .. code-block:: python

        class MyOption(Option):
            ""Docstring for MyOption""

            A = "a", "Option A"
            B = "b", "Option B"
            C = "c", "Option C"

            @staticmethod
            def create(option: str):
                return Option._create(MyOption, option)

            @staticmethod
            def options(include_docs: bool = False):
                return Option._options(MyOption, include_docs=include_docs)

    The resulting class can then be used as follows:

    .. code-block:: python

        assert MyOption.A == "a"

        assert MyOption.create("a") == MyOption.A
        assert MyOption.create("a") == "a"

        assert MyOption.options() == ["a", "b", "c"]
        assert MyOption.options(include_docs=True) == [
            ("a", "Option A"),
            ("b", "Option B"),
            ("c", "Option C"),
        ]
    """

    def __new__(cls, value, doc=None):
        self = super().__new__(cls, value)
        self._value_ = value

        if doc is not None:
            self.__doc__ = doc

        return self

    def __str__(self) -> str:
        return str(self.value)

    def __eq__(self, other) -> bool:
        return str(self) == str(other)

    @staticmethod
    def canonicalise(option: str) -> str:
        """
        Return the canonical version of the passed option string.
        This is a lower case string with no extra whitespace (none
        at beginning or end, with repeated whitespace removed),
        and where hyphens are replaced by underscores
        """
        option = str(option)
        option = option.replace("-", "_")
        option = option.lower().lstrip().rstrip()
        option = " ".join(option.split())
        return option

    @staticmethod
    def _create(CLS, option):
        """
        Create an option from the passed string. Returns
        an option of type CLS
        """
        if type(option) is CLS:
            return option
        else:
            # Clean the option up, lowercasing, removing whitespace etc.
            option = CLS.canonicalise(option)

            try:
                return CLS(option)
            except ValueError:
                raise ValueError(
                    f"Invalid option '{option}'. Available options are: "
                    + ", ".join(Option._options(CLS))
                )

    @staticmethod
    def _options(CLS, include_docs: bool = False):
        """
        Return the list of strings representing all
        of the available options
        """
        if include_docs:
            return [(x.value, x.__doc__) for x in CLS]
        else:
            return [x.value for x in CLS]
