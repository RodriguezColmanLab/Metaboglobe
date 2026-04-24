from metaboglobe.plotting._vector_2d import Vector2


class Box2:
    """Represents a 2D bounding box. The object is immutable, and can have a width/height of zero. (Both not negative.)
    """
    _min: Vector2
    _max: Vector2

    @staticmethod
    def enclosing(*args: Vector2) -> "Box2":
        """Gets a Box2 that encloses all given points. So `box.max.x == max(x of args)`, and same for `min` and `y`."""
        min_pos = Vector2.min_x_y(*args)
        max_pos = Vector2.max_x_y(*args)
        return Box2(min_pos, max_pos)

    def __init__(self, min: Vector2, max: Vector2):
        """Initializes a Box2 with the given bounding box. Raises a ValueError if min.x > max.x or if min.y > max.y."""
        if min.x > max.x:
            raise ValueError("min x must be less than max")
        if min.y > max.y:
            raise ValueError("min y must be less than max")
        self._min = min
        self._max = max

    @property
    def min(self) -> Vector2:
        """The top left corner point of the box."""
        return self._min

    @property
    def max(self) -> Vector2:
        """The bottom right corner point of the box. Both the x and y are greater than or equal to the min."""
        return self._max

    def width(self) -> float:
        """The width of the box, zero or larger."""
        return self._max.x - self._min.x

    def height(self) -> float:
        """The height of the box, zero or larger."""
        return self._max.y - self._min.y

    def area(self) -> float:
        """The area of the box, zero or larger."""
        return self.width() * self.height()





