import math
from typing import Literal, NamedTuple

import numpy
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.cbook import simple_linear_interpolation
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle, PathPatch
from matplotlib.path import Path
from matplotlib.text import Text
from matplotlib.transforms import Bbox
from numpy import ndarray

from metaboglobe.plotting._text_size import estimate_width_height


class TextWithAnchor(NamedTuple):
    """Holds the text to draw."""
    x: float
    y: float
    text: str


class _TextBbox(NamedTuple):
    """Holds an estimation of the size of text."""
    x: float
    y: float
    width: float
    height: float
    ha: Literal["center", "left", "right"]
    va: Literal["center", "top", "bottom"]

    @staticmethod
    def from_artist(text: Text) -> "_TextBbox":
        x, y = text.get_position()
        width, height = estimate_width_height(text.get_text(), text.get_fontsize())
        return _TextBbox(
            x=x, y=y,
            width=width, height=height,
            ha=text.get_horizontalalignment(), va=text.get_verticalalignment()
        )


    def to_bbox(self) -> Bbox:
        if self.ha == "left":
            x_start = self.x
        elif self.ha == "right":
            x_start = self.x - self.width
        else:
            x_start = self.x - self.width / 2

        if self.va == "top":
            y_start = self.y
        elif self.va == "bottom":
            y_start = self.y - self.height
        else:
            y_start = self.y - self.height / 2
        return Bbox.from_bounds(x_start, y_start, self.width, self.height)


class CollisionMap:

    _x_offset: int
    _y_offset: int
    _resolution: float
    _grid: ndarray  # On a scale of 0 to 255

    def __init__(self, ax: Axes, resolution: float = 5):
        """Initializes an empty collision map that spans the given axis."""
        self._resolution = resolution

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        self._x_offset = math.floor(min(xlim))
        self._y_offset = math.floor(min(ylim))
        max_x = max(xlim)
        max_y = max(ylim)
        grid_width = int(math.ceil((max_x - self._x_offset) / resolution))
        grid_height = int(math.ceil((max_y - self._y_offset) / resolution))
        self._grid = numpy.zeros((grid_height, grid_width), dtype=numpy.uint8)

    def _mark_on_grid(self, x: float, y: float, width: float = 0, height: float = 0, weight: float = 1.0):
        grid_x = int((x - self._x_offset) / self._resolution)  # Inclusive
        grid_y = int((y - self._y_offset) / self._resolution)  # Inclusive

        grid_max_x = int((x + width - self._x_offset) / self._resolution) + 1  # Exclusive
        grid_max_y = int((y + height - self._y_offset) / self._resolution) + 1  # Exclusive

        if grid_max_x <= 0 or grid_x >= self._grid.shape[1] or grid_max_y <= 0 or grid_y >= self._grid.shape[0]:
            return

        # Clip to boundaries of grid
        grid_x = max(grid_x, 0)
        grid_y = max(grid_y, 0)
        grid_max_x = min(grid_max_x, self._grid.shape[1])
        grid_max_y = min(grid_max_y, self._grid.shape[0])

        self._grid[grid_y:grid_max_y, grid_x:grid_max_x] = int(weight * 255)

    def _fraction_free(self, location: _TextBbox) -> float:
        """Checks if the given position is free on the grid. (True if free, False if occupied.)"""
        bbox = location.to_bbox()

        # All values are inclusive here
        grid_x_min = int((bbox.xmin - self._x_offset) / self._resolution)
        grid_x_max = int((bbox.xmax - self._x_offset) / self._resolution)
        grid_y_min = int((bbox.ymin - self._y_offset) / self._resolution)
        grid_y_max = int((bbox.ymax - self._y_offset) / self._resolution)

        total_count = 0
        occupied_count = 0
        for grid_x in range(grid_x_min, grid_x_max + 1):
            for grid_y in range(grid_y_min, grid_y_max + 1):
                if grid_x < 0 or grid_y < 0 or grid_x > self._grid.shape[1] or grid_y > self._grid.shape[0]:
                    occupied_count += 1  # Out of bounds, consider as occupied
                occupied_count += self._grid[grid_y, grid_x] / 255  # Occupied
                total_count += 1
        return (total_count - occupied_count) / total_count

    def add_artist(self, artist: Artist):
        """Adds an artist to the collision map. We're selective in which artists we support, since supporting all of
        them would be too much work."""
        if isinstance(artist, Text):
            # We just mark the center position
            if not artist.get_text() or artist == artist.axes.title:
                return  # Ignore empty strings
            bbox = _TextBbox.from_artist(artist).to_bbox()

            self._mark_on_grid(bbox.xmin, bbox.ymin, bbox.width, bbox.height)
        elif isinstance(artist, FancyArrowPatch) or isinstance(artist, PathPatch):
            # Plotting over reaction arrows is a bigger problem than over pathway links
            weight = 1.0 if isinstance(artist, FancyArrowPatch) else 0.5

            # We have to follow the path
            path = artist.get_path()
            vertices, codes = path.vertices, path.codes
            if codes is None:  # If no codes have been specified, generate them ourselves
                codes = [Path.LINETO] * len(vertices)
                codes[0] = Path.MOVETO

            previous_vertex = None
            for vertex, code in zip(vertices, codes):
                if code == Path.MOVETO or code == Path.CLOSEPOLY:
                    previous_vertex = vertex
                    continue  # Not supported
                if previous_vertex is not None:
                    length = math.sqrt((vertex[0] - previous_vertex[0]) ** 2 + (vertex[1] - previous_vertex[1]) ** 2)
                    steps = math.ceil(length / self._resolution)
                    if steps > 1:
                        interpolated = simple_linear_interpolation(numpy.array([vertex, previous_vertex]), steps)
                    else:
                        interpolated = [vertex, previous_vertex]
                    for interpolated_vertex in interpolated:
                        self._mark_on_grid(*interpolated_vertex, weight=weight)
                previous_vertex = vertex
        elif isinstance(artist, Circle):
            # We just mark the center position
            x, y = artist.get_center()
            radius_reduced = artist.get_radius() * 0.8
            self._mark_on_grid(x - radius_reduced, y - radius_reduced, radius_reduced * 2, radius_reduced * 2)
        elif isinstance(artist, Rectangle):
            # We mark the entire rectangle
            self._mark_on_grid(artist.get_x(), artist.get_y(), width=artist.get_width(), height=artist.get_height())


    def fit_text(self, ax: Axes, texts: list[TextWithAnchor], *, fontsize: float = 8, zorder: int = 0):
        # Sort to have longest texts first (theses are the most difficult to find a spot for)
        texts.sort(key=lambda t: len(t.text), reverse=True)
        for text_with_anchor in texts:
            text = text_with_anchor.text
            x = text_with_anchor.x
            y = text_with_anchor.y

            text_width, text_height = estimate_width_height(text, fontsize)

            highest_score = 0
            highest_scoring_option = None
            for option, desirability in {
                    _TextBbox(x - 10, y, text_width, text_height, "right", "center"): 1.0,
                    _TextBbox(x + 10, y, text_width, text_height, "left", "center"): 0.99,
                    _TextBbox(x, y - 17, text_width, text_height, "center", "top"): 0.99,
                    _TextBbox(x, y + 17, text_width, text_height, "center", "bottom"): 0.99,
                    _TextBbox(x - 6, y - 17, text_width, text_height, "right", "top"): 0.98,
                    _TextBbox(x + 6, y + 19, text_width, text_height, "left", "bottom"): 0.98,
                    _TextBbox(x + 6, y - 17, text_width, text_height, "left", "top"): 0.98,
                    _TextBbox(x - 6, y + 19, text_width, text_height, "right", "bottom"): 0.98,
                    _TextBbox(x, y, text_width, text_height, "center", "center"): 0.5  # Last resort
                }.items():
                score = self._fraction_free(option) + desirability
                if score > highest_score:
                    highest_score = score
                    highest_scoring_option = option

            # Draw the best option
            option = highest_scoring_option
            ax.text(option.x, option.y, text, fontsize=fontsize, ha=option.ha, va=option.va, zorder=zorder)

            # Block on the map
            bbox = option.to_bbox()
            self._mark_on_grid(bbox.xmin, bbox.ymin, bbox.width, bbox.height)


    def debug_draw(self, ax: Axes):
        ax.imshow(self._grid, extent=(self._x_offset, self._x_offset + self._resolution * self._grid.shape[1],
                                      self._y_offset + self._resolution * self._grid.shape[0], self._y_offset),
                  cmap="spring", vmin=0, vmax=255, zorder=-10)
