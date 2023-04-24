__all__ = ["ProgressBar"]

from ..legacy.Base import ProgressBar as _ProgressBar


def _in_notebook():
    """
    https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    """
    try:
        from IPython import get_ipython

        if "IPKernelApp" not in get_ipython().config:  # pragma: no cover
            return False
    except ImportError:
        return False
    except AttributeError:
        return False
    return True


_forced_jupyter = False


def _force_jupyter():
    global _forced_jupyter

    if _forced_jupyter:
        return

    if _in_notebook():
        _ProgressBar.set_theme("jupyter")

    _forced_jupyter = True


class ProgressBar:
    def __init__(self, total=None, text=None, _bar=None):
        _force_jupyter()

        if _bar is not None:
            self._bar = _bar
        elif total is not None:
            if text is not None:
                self._bar = _ProgressBar(total=total, text=text)
            else:
                self._bar = _ProgressBar(total=total)
        else:
            if text is not None:
                self._bar = _ProgressBar(text=text)
            else:
                self._bar = _ProgressBar()

    def __enter__(self):
        self._bar = self._bar.enter()
        return self

    def __exit__(self, type, value, tb):
        self._bar.exit()

    def tick(self, text=None):
        if text is not None:
            self._bar.tick(text)
        else:
            self._bar.tick()

    def set_speed_unit(self, unit):
        self._bar.set_speed_unit(unit)

    def set_progress(self, progress, text=None):
        if text is not None:
            self._bar.set_progress(progress, text)
        else:
            self._bar.set_progress(progress)

    def success(self, text=None):
        if text is None:
            self._bar.success()
        else:
            self._bar.success(text)

    def failure(self, text=None):
        if text is None:
            self._bar.failure()
        else:
            self._bar.failure(text)

    @staticmethod
    def set_theme(theme):
        if (not _in_notebook()) or (theme.lower().find("silent") != -1):
            _ProgressBar.set_theme(theme)

    @staticmethod
    def set_silent():
        _ProgressBar.set_silent()
