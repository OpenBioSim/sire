from dataclasses import dataclass as _dataclass
from dataclasses import field as _field
from typing import List as _List
from typing import Dict as _Dict

from ._console import Console

__all__ = ["Simple"]


@_dataclass
class Simple:
    """This is the monochrome 'Simple' theme"""

    frames: _Dict[int, _List[str]] = _field(default_factory=dict)

    def should_highlight(self):
        return False

    def highlighter(self):
        return None

    def should_markup(self):
        return False

    def panel_box(self, style):
        from rich import box as _box

        return _box.SQUARE

    def padding_style(self, style):
        return ""

    def text(self, style):
        return "white"

    def error(self):
        return "white"

    def warning(self):
        return "white"

    def info(self):
        return "white"

    def error_text(self):
        return f"bold {self.error()}"

    def warning_text(self):
        return f"bold {self.warning()}"

    def info_text(self):
        return f"bold {self.info()}"

    def spinner(self, text: str):
        console = Console._get_console()
        return console.status(text, spinner="dots")

    def rule(self, style):
        return "white"

    def panel(self, style):
        return "white"
