import numpy
from matplotlib.path import Path

from metaboglobe.plotting._vector_2d import Vector2


class Curve2:
    """A mutable curve."""

    _vertices: list[Vector2]
    _codes: list[numpy.unsignedinteger]

    def __init__(self, start: Vector2):
        self._vertices = [start]
        self._codes = [Path.MOVETO]


    def to_path(self) -> Path:
        """Converts the curve to a Path."""
        return Path(self._vertices, self._codes)

    def append_line_to(self, point: Vector2) -> None:
        """Appends a line to the given point."""
        self._vertices.append(point)
        self._codes.append(Path.LINETO)

    def append_rounded_corner_to(self, *, corner: Vector2, end: Vector2, radius_px: float = 8) -> None:
        """Draws a line to the given end point, with a rounded corner of the given radius around the given corner
        location. So you'll get a straight line, until we come into the radius around the corner, where the curve
        starts. Once outside the radius around the corner, the line resumes as a straight line towards the end.

        The radius is automatically shrunken if it's larger than the distance from the current position to the corner,
        or from the corner to the end.
        """
        start = self._vertices[-1]

        # Determine the maximum possible radius
        radius_px = min(start.distance(corner) * 0.99, corner.distance(end) * 0.99, radius_px)
        if radius_px < 0:
            raise ValueError("radius_px cannot be negative")

        corner_start_point = corner.towards(start, radius_px)
        corner_end_point = corner.towards(end, radius_px)

        self._vertices.append(corner_start_point)
        self._codes.append(Path.LINETO)

        self._vertices.append(corner)
        self._codes.append(Path.CURVE3)

        self._vertices.append(corner_end_point)
        self._codes.append(Path.CURVE3)

        self._vertices.append(end)
        self._codes.append(Path.LINETO)

    def append_line_then_curve_to(self, goal: Vector2, *, curved_distance: float = 20) -> None:
        """Appends a straight line to the given goal point. At curved_distance from the goal, the straight line ends
        and becomes a curve instead."""
        start = self._vertices[-1]

        curved_distance = min(curved_distance, start.distance(goal) * 0.99)  # Length of the rounded part of the arrow
        self._vertices.append(goal.towards(start, curved_distance))
        self._codes.append(Path.LINETO)

        self._vertices.append(goal)
        self._codes.append(Path.CURVE3)

    def append_curve_then_line_to(self, goal: Vector2, *, curved_distance: float = 20) -> None:
        """Goes into a curve towards the given goal for curved_distance. From then on, we draw a straight line to the
        goal."""
        start = self._vertices[-1]

        curved_distance = min(curved_distance, start.distance(goal) * 0.99)  # Length of the rounded part of the arrow
        self._vertices.append(start.towards(goal, curved_distance))
        self._codes.append(Path.CURVE3)
        self._vertices.append(goal)
        self._codes.append(Path.LINETO)

    def shorten(self, length_start: float, length_end: float) -> None:
        """Removes the given length from the start and end of the curve."""
        while length_start > 0 and len(self._vertices) >= 2:
            current_segment_length = self._vertices[0].distance(self._vertices[1])
            if current_segment_length > length_start:
                # Remove part of segment
                new_start = self._vertices[1].towards(self._vertices[0], length_start)
                self._vertices[0] = new_start
                length_start = 0
            else:
                # Remove entire segment
                length_start -= current_segment_length
                del self._vertices[0]
                del self._codes[0]
                self._codes[0] = Path.MOVETO

        while length_end > 0 and len(self._vertices) >= 2:
            current_segment_length = self._vertices[0].distance(self._vertices[1])
            if current_segment_length > length_end:
                # Remove part of segment
                new_end = self._vertices[-2].towards(self._vertices[-1], length_end)
                self._vertices[-1] = new_end
                length_end = 0
            else:
                # Remove entire segment
                length_end -= current_segment_length
                del self._vertices[-1]
                del self._codes[-1]