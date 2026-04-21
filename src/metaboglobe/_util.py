import itertools
import math
from typing import Iterable

_ENANTIOMER_MARKERS = ["alpha-", "beta-", "gamma-", "D-", "L-"]


MPLColor = tuple[float, float, float] \
           | str \
           | tuple[float, float, float, float] \
           | tuple[tuple[float, float, float] | str, float]


class Direction:
    """Represents a unit vector in 2D in a given direction."""

    __slots__ = "_angle_radians"

    _angle_radians: float

    def __init__(self, angle_radians: float):
        self._angle_radians = angle_radians % (2 * math.pi)

    @property
    def angle_radians(self) -> float:
        return self._angle_radians

    @property
    def angle_degrees(self) -> float:
        return math.degrees(self._angle_radians)

    def orthogonal(self) -> "Direction":
        """Returns a direction orthogonal to this one. (90 degrees to the right.)"""
        return Direction(angle_radians=self.angle_radians + math.pi / 2)

    def opposite(self) -> "Direction":
        """Returns a direction opposite to this one."""
        return Direction(angle_radians=self.angle_radians + math.pi)

    def dx(self) -> float:
        """Returns the x component of the direction."""
        return math.cos(self._angle_radians)

    def dy(self) -> float:
        """Returns the y component of the direction."""
        return math.sin(self._angle_radians)

    def __repr__(self) -> str:
        return f"Direction({self.angle_degrees:.1f})"


def point_direction(x1: float, y1: float, x2: float, y2: float) -> Direction:
    """Gets the direction from point 1 to point 2."""
    angle = math.atan2(y2 - y1, x2 - x1)
    return Direction(angle_radians=angle)



def optimize_for_display(name: str) -> str:
    """Converts a KEGG compound name like "alpha-D-glucose-6-phosphate" to a display name suitable for Matplotlib
    like "$\\alpha$-D-glucose-6P"."""
    name = name.replace("-phosphate", "P")
    name = name.replace("-bisphosphate", "P$_2$")
    name = name.replace("alpha", "$\\alpha$")
    name = name.replace("beta", "$\\beta$")
    name = name.replace(" ", "-")
    return name


def optimize_for_matching(name: str) -> str:
    return name.lower().replace(" ", "-")


def wrap_text(input_str: str, max_line_length: int) -> str:
    """Wraps the input string to a maximum line length, breaking at spaces."""
    words = input_str.split()
    lines = []
    current_line = ""

    for word in words:
        if len(current_line) + len(word) + 1 <= max_line_length:
            if current_line:
                current_line += " "
            current_line += word
        else:
            lines.append(current_line)
            current_line = word

    if current_line:
        lines.append(current_line)

    return "\n".join(lines)


def _stereoisomer_search_at_start(haystack: str) -> int:
    """Returns 0 if no enantiomer start (such as "D-") is found. Otherwise, returns the length of the enantiomer start.
    """
    for enantiomer_marker in _ENANTIOMER_MARKERS:
        if haystack.startswith(enantiomer_marker):
            return len(enantiomer_marker)
    return 0


def get_names_without_stereoisomers(name: str) -> Iterable[str]:
    """Removes enantiomer markers like "D-" and "alpha-" from the name. Returns all possible variations without at least
    one specifier. So "alpha-D-glucose-6-phosphate" would return ["alpha-glucose-6-phosphate", "D-glucose-6-phosphate",
    "glucose-6-phosphate"]."""
    stereoisomers_all = set()

    i = 0
    while i < len(name):
        if i == 0 or name[i - 1] in "-(":
            # Start of a new word, check for enantiomer marker
            name_from_here = name[i:]
            stereoisomer_length = _stereoisomer_search_at_start(name_from_here)
            if stereoisomer_length > 0:
                stereoisomers_all.add((i, i + stereoisomer_length))  # Found an stereoisomer
                i += stereoisomer_length
                continue

        i += 1

    # Find all possible subsets of stereoisomer_variant, excluding the empty set
    for stereroisomer_variant in _all_possible_subsets(stereoisomers_all):
        name_without_enantiomers = name
        for start, end in sorted(stereroisomer_variant, reverse=True):
            name_without_enantiomers = name_without_enantiomers[:start] + name_without_enantiomers[end:]
        yield name_without_enantiomers


def _all_possible_subsets(s: set) -> list[set]:
    """Returns a list of all possible subsets of the given set, excluding the empty set."""
    s = list(s)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(1, len(s)+1))
