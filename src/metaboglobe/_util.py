import itertools
import math
from typing import Iterable

_ENANTIOMER_MARKERS = ["alpha-", "beta-", "gamma-", "D-", "L-"]


MPLColor = tuple[float, float, float] \
           | str \
           | tuple[float, float, float, float] \
           | tuple[tuple[float, float, float] | str, float]





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
