
_CHARACTER_WIDTHS = {
    # Most characters are assumed to have width of 1. Some longer sequences also have 1, which are listed in this map,
    # and some characters are smaller, which are also listed in this map.
    "$\\alpha$": 1,
    "$\\beta$": 1,
    "$\\gamma$": 1,
    "$\\delta$": 1,
    "$_2$": 0.6,
    ".": 0.5,
    ",": 0.5,
    "(": 0.5,
    ")": 0.5,
}


def estimate_width_height(text: str, fontsize: float = 8) -> tuple[float, float]:
    """Estimates the width and height of the given text, based on the font size. It's a simple calculation, but does
    know that for example $\\alpha$ and "a" have the same width."""
    lines = text.split("\n")
    text_height = len(lines) * fontsize * 1.2
    text_width = max((_estimate_line_width(line, fontsize) for line in lines), default=0)

    return text_width, text_height


def _estimate_line_width(text: str, fontsize: float = 8) -> float:
    i = 0
    total_width = 0
    while i < len(text):
        for replacement, width in _CHARACTER_WIDTHS.items():
            replacement_length = len(replacement)
            if i + replacement_length > len(text):
                continue
            if text[i:i + replacement_length] == replacement:
                total_width += width
                i += replacement_length
                break

        # No match found, just normal increase
        total_width += 1
        i += 1

    return total_width * fontsize * 0.9
