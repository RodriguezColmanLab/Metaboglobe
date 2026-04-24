import numpy
from matplotlib.axes import Axes
from matplotlib.path import Path

from metaboglobe.plotting._vector_2d import Vector2, Direction


def _average(direction1: Direction | None, direction2: Direction | None) -> Direction:
    if direction1 is None:
        if direction2 is None:
            raise ValueError("direction1 and direction2 cannot both be None")
        return direction2
    if direction2 is None:
        return direction1

    return direction1.middle(direction2)


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

    def append_cut_corner_to(self, corner: Vector2, end: Vector2) -> None:
        """Draws a line to the given end point, cutting out a diagonal corner as large as possible."""
        start = self._vertices[-1]
        corner_size = min(start.distance(corner) * 0.99, corner.distance(end) * 0.99)

        corner_start = corner.towards(start, corner_size)
        corner_end = corner.towards(end, corner_size)

        self._vertices.append(corner_start)
        self._codes.append(Path.LINETO)
        self._vertices.append(corner_end)
        self._codes.append(Path.LINETO)
        self._vertices.append(end)
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

    def split(self, separation_distance: float = 4) -> tuple["Curve2", "Curve2"]:
        """For reversible reactions, this duplicates the curve into two, one for the forwards reaction, one for the
        backwards. Both lines run parallel to each other with some space in between, in opposite direction."""
        vertices_a = list()
        vertices_b = list()

        for i, vertex in enumerate(self._vertices):
            direction_forwards = vertex.direction(self._vertices[i + 1]) if i < len(self._vertices) - 1 else None
            direction_backwards = vertex.direction(self._vertices[i - 1]).opposite() if i > 0 else None
            direction = _average(direction_forwards, direction_backwards)

            offset = direction.orthogonal().to_vector(separation_distance / 2)
            vertices_a.append(vertex + offset)
            vertices_b.append(vertex - offset)

        vertices_b.reverse()

        path_a = Curve2(vertices_a[0])
        path_a._vertices = vertices_a
        path_a._codes = self._codes.copy()

        codes = self._codes.copy()
        del codes[0]
        codes.append(Path.MOVETO)
        codes.reverse()

        path_b = Curve2(vertices_b[0])
        path_b._vertices = vertices_b
        path_b._codes = codes

        return path_a, path_b

    def shorten_both_sides(self, shorten_distance: float):
        """Shortens the curve with the given distance on both sides."""
        self._shorten_front(shorten_distance)
        self._shorten_back(shorten_distance)

    def _shorten_front(self, shorten_distance: float):
        while shorten_distance > 0 and len(self._vertices) >= 2:
            length_first_segment = self._vertices[0].distance(self._vertices[1])

            # Remove this segment entirely
            if length_first_segment <= shorten_distance:
                del self._vertices[0]
                del self._codes[0]
                self._codes[0] = Path.MOVETO

                shorten_distance -= length_first_segment
                continue

            # Need to shorten a line
            if length_first_segment > shorten_distance:
                self._vertices[0] = self._vertices[1].towards(self._vertices[0], length_first_segment - shorten_distance)

                if len(self._vertices) >= 4 and self._codes[-1] == Path.CURVE3 and self._codes[-2] == Path.CURVE3:
                    # Need to adapt the curve, otherwise the arrow pointer looks strange

                    # Imagine the situation like this:
                    # [0] is MOVETO, [1] and [2] are CURVE3, and [3] is LINETO
                    #
                    # [1]    _[2]-------------------------------------------------> [3]
                    #      /
                    # [0]|
                    #
                    # We do a similar operation as for _shorten_back: we move point [1] a lot closer to point [2]

                    distance_1_2 = self._vertices[1].distance(self._vertices[2])
                    new_1 = self._vertices[2].towards(self._vertices[1], distance_1_2 * 0.33)
                    self._vertices[1] = new_1

                return  # We're done shortening

            return  # Ran out of options to shorten

    def extend_both_sides(self, extend_distance: float):
        """Extends the curve with the given distance on both sides, by either extendind the straight line at the end,
        or by appending a straight line."""

        if len(self._vertices) < 2:
            return

        # Extend at start
        if self._codes[1] == Path.LINETO:
            # Can simply extend this line
            self._vertices[0] = self._vertices[0].towards(self._vertices[1], -extend_distance)
        else:
            # Need to append a line segment
            self._vertices.insert(0, self._vertices[0].towards(self._vertices[1], -extend_distance))
            self._codes.insert(1, Path.LINETO)  # Index 0 is always Path.MOVETO, so use index 1

        # Extend at end
        if self._codes[-1] == Path.LINETO:
            # Can simply extend this line
            self._vertices[-1] = self._vertices[-1].towards(self._vertices[-2], -extend_distance)
        else:
            # Need to append a line segment
            self._vertices.append(self._vertices[-1].towards(self._vertices[-2], -extend_distance))
            self._codes.append(Path.LINETO)


    def _shorten_back(self, shorten_distance: float):
        while shorten_distance > 0 and len(self._vertices) >= 2:
            length_last_segment = self._vertices[-1].distance(self._vertices[-2])

            # Remove this segment entirely
            if length_last_segment <= shorten_distance:
                del self._vertices[-1]
                del self._codes[-1]

                shorten_distance -= length_last_segment
                continue

            # Need to shorten a line
            if length_last_segment > shorten_distance:
                self._vertices[-1] = self._vertices[-2].towards(self._vertices[-1], length_last_segment - shorten_distance)

                if len(self._vertices) >= 4 and self._codes[-1] == Path.CURVE3 and self._codes[-2] == Path.CURVE3:
                    # Need to adapt the curve, otherwise the arrow pointer looks strange

                    # Imagine the situation like this:
                    # [-1] and [-2] are CURVE3, [-3] and [-4] are LINETO
                    #
                    # [-2]    _[-3]------------------------------------------------- [-4]
                    #      /
                    # [-1]V
                    #
                    # Here, the arrow points down at point [-1], because of point [-2]. With the shortened line between
                    # [-2] and [-1], that looks unfortunate, as the arrow head is longer than the line [-2][-1].
                    # So we move [-2] closer to [-3]

                    distance_2_3 = self._vertices[-2].distance(self._vertices[-3])
                    new_2 = self._vertices[-3].towards(self._vertices[-2], distance_2_3 * 0.33)
                    self._vertices[-2] = new_2

                return  # We're done shortening

            return  # Ran out of options to shorten

    def debug_draw(self, ax: Axes):
        """Draws the vertex points of the curve, colored by their type (start, line, curve). This is only for debugging
         purposes."""
        x_positions = [vertex.x for vertex in self._vertices]
        y_positions = [vertex.y for vertex in self._vertices]

        colors = list()
        for code in self._codes:
            if code == Path.MOVETO or code == Path.LINETO:
                colors.append(0)
            else:
                colors.append(1)

        ax.scatter(x_positions, y_positions, c=colors, s=16, lw=0, cmap="winter")
