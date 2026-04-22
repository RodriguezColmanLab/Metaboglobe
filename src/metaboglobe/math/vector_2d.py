import math
from typing import NamedTuple, Iterable

from numpy import sign



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


class Vector2(NamedTuple):
    """Represents a vector in 2D space."""

    x: float
    y: float

    @staticmethod
    def min_x_y(*args: "Vector2") -> "Vector2":
        """Returns a vector with the minimum x and y of the given vectors."""
        min_x = args[0].x
        min_y = args[0].y
        for other in args:
            min_x = min(min_x, other.x)
            min_y = min(min_y, other.y)
        return Vector2(x=min_x, y=min_y)

    @staticmethod
    def max_x_y(*args: "Vector2") -> "Vector2":
        """Returns a vector with the maximum x and y of the given vectors."""
        max_x = args[0].x
        max_y = args[0].y
        for other in args:
            max_x = max(max_x, other.x)
            max_y = max(max_y, other.y)
        return Vector2(x=max_x, y=max_y)

    def with_x(self, x: float) -> "Vector2":
        """Returns a copy of this vector with the given x component."""
        if x == self.x:
            return self
        return Vector2(x=x, y=self.y)

    def with_y(self, y: float) -> "Vector2":
        """Returns a copy of this vector with the given y component."""
        if y == self.y:
            return self
        return Vector2(x=self.x, y=y)

    def towards(self, goal: "Vector2", distance: float) -> "Vector2":
        """Returns a copy of this vector that is the given distance away from this point, into the direction of the
        goal. Note that if the goal is closer than the given distance, this method will overshoot."""
        dx = goal.x - self.x
        dy = goal.y - self.y
        if dx == 0:
            return Vector2(x=self.x, y=self.y + distance * sign(dy))
        elif dy == 0:
            return Vector2(x=self.x + distance * sign(dx), y=self.y)
        else:
             # More complicated math
            direction = self.direction(goal)
            return Vector2(x=self.x + direction.dx() * distance, y=self.y + direction.dy() * distance)

    def direction(self, other: "Vector2") -> Direction:
        """Gets the direction from this point to the given point."""
        return point_direction(self.x, self.y, other.x, other.y)

    def distance(self, other: "Vector2") -> float:
        """Gets the distance to the other point."""
        dx = self.x - other.x
        dy = self.y - other.y
        if dx == 0:
            return abs(dy)
        elif dy == 0:
            return abs(dx)
        else:
            return math.sqrt(dx * dx + dy * dy)


def point_direction(x1: float, y1: float, x2: float, y2: float) -> Direction:
    """Gets the direction from point 1 to point 2."""
    angle = math.atan2(y2 - y1, x2 - x1)
    return Direction(angle_radians=angle)