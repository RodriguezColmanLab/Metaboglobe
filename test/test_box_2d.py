from unittest import TestCase

from metaboglobe.math.box_2d import Box2
from metaboglobe.math.vector_2d import Vector2


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

