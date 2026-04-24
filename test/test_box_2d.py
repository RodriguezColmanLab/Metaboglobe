from unittest import TestCase

from metaboglobe.plotting._box_2d import Box2
from metaboglobe.plotting._vector_2d import Vector2


class Test(TestCase):

    def test_creation(self):
        corner_a = Vector2(10, 10)
        corner_b = Vector2(12, 11)
        box = Box2(corner_a, corner_b)
        self.assertEqual(box.min, corner_a)
        self.assertEqual(box.max, corner_b)

    def test_creation_rejects_inverted_x_range(self):
        with self.assertRaisesRegex(ValueError, "min x must be less than max"):
            Box2(min=Vector2(2, 0), max=Vector2(1, 1))

    def test_creation_rejects_inverted_y_range(self):
        with self.assertRaisesRegex(ValueError, "min y must be less than max"):
            Box2(min=Vector2(0, 2), max=Vector2(1, 1))

    def test_creation_allows_flat_box(self):
        box = Box2(min=Vector2(3, 3), max=Vector2(3, 3))
        self.assertEqual(box.width(), 0)
        self.assertEqual(box.height(), 0)
        self.assertEqual(box.area(), 0)

    def test_enclosing(self):
        box = Box2.enclosing(
            Vector2(3, 4),
            Vector2(-1, 6),
            Vector2(2, -2),
            Vector2(5, 5),
        )
        self.assertEqual(box.min, Vector2(-1, -2))
        self.assertEqual(box.max, Vector2(5, 6))

    def test_width_height_area(self):
        box = Box2(min=Vector2(1, 2), max=Vector2(6, 7))
        self.assertEqual(box.width(), 5)
        self.assertEqual(box.height(), 5)
        self.assertEqual(box.area(), 25)

    def test_is_on_border_true_for_edges_and_corners(self):
        box = Box2(min=Vector2(1, 1), max=Vector2(4, 3))
        border_points = [
            Vector2(1, 1),  # bottom-left corner
            Vector2(4, 1),  # bottom-right corner
            Vector2(1, 3),  # top-left corner
            Vector2(4, 3),  # top-right corner
            Vector2(2, 1),  # bottom edge
            Vector2(2, 3),  # top edge
            Vector2(1, 2),  # left edge
            Vector2(4, 2),  # right edge
        ]
        for point in border_points:
            self.assertTrue(box.is_on_border(point), msg=f"Expected border point: {point}")

    def test_is_on_border_false_for_inside_and_outside_points(self):
        box = Box2(min=Vector2(1, 1), max=Vector2(4, 3))
        non_border_points = [
            Vector2(2, 2),  # inside
            Vector2(0, 2),  # outside left
            Vector2(5, 2),  # outside right
            Vector2(2, 0),  # outside below
            Vector2(2, 4),  # outside above
            Vector2(0, 0),  # outside diagonally
        ]
        for point in non_border_points:
            self.assertFalse(box.is_on_border(point), msg=f"Expected non-border point: {point}")
