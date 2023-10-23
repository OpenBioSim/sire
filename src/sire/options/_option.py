__all__ = ["Option"]

try:
    from enum import StrEnum as _StrEnum
except ImportError:
    # define a basic StrEnum class if it is not available
    from enum import Enum

    class _StrEnum(str, Enum):
        def __str__(self) -> str:
            return str(self.value)


class Option(_StrEnum):
    """
    Base class of all of the Option objects. These are
    effectively StrEnum objects that have additional
    functions to support case-insentive assigment
    and querying
    """

    @staticmethod
    def canonicalise(option: str) -> str:
        """
        Return the canonical version of the passed option string
        """
        option = str(option)
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
    def _options(CLS):
        """
        Return the list of strings representing all
        of the available options
        """
        return [x.value for x in CLS]
