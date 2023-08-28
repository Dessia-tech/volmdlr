"""
Showcase of performing boolean operations using pixelization.
"""

from volmdlr import Point2D, Vector2D
from volmdlr.discrete_representation import PointBasedPixelization
from volmdlr.wires import ClosedPolygon2D

PIXEL_SIZE = 0.05

# Define a polygon
points = [
    Point2D(5, 1),
    Point2D(5.25, 0.5),
    Point2D(6, 0.5),
    Point2D(5.45, 0),
    Point2D(6, -1),
    Point2D(5, -0.5),
    Point2D(4, -1),
    Point2D(4.55, 0),
    Point2D(4, 0.5),
    Point2D(4.75, 0.5),
    Point2D(5, 1),
]

closed_polygon_1 = ClosedPolygon2D(points)
closed_polygon_2 = closed_polygon_1.translation(Vector2D(1.2, -0.5))

pixelization_1 = PointBasedPixelization.from_closed_polygon(closed_polygon_1, PIXEL_SIZE).fill_enclosed_pixels()
pixelization_2 = PointBasedPixelization.from_closed_polygon(closed_polygon_2, PIXEL_SIZE).fill_enclosed_pixels()

# Union
union = pixelization_1 + pixelization_2  # equivalent to pixelization1.union(pixelization_2)
ax = union.plot(color="blue")
closed_polygon_1.plot(ax=ax)
closed_polygon_2.plot(ax=ax)

# Difference
difference_1 = pixelization_1 - pixelization_2  # equivalent to pixelization1.difference(pixelization_2)
difference_2 = pixelization_2 - pixelization_1  # equivalent to pixelization2.difference(pixelization_1)
ax = difference_1.plot(color="purple")
difference_2.plot(ax=ax, color="red")
closed_polygon_1.plot(ax=ax)
closed_polygon_2.plot(ax=ax)

# Intersection
intersection = pixelization_1 and pixelization_2  # equivalent to pixelization1.intersection(pixelization_2)
ax = intersection.plot(color="green")
closed_polygon_1.plot(ax=ax)
closed_polygon_2.plot(ax=ax)

# Symmetric difference
symmetric_diff = pixelization_1 ^ pixelization_2  # equivalent to pixelization1.symmetric_difference(pixelization_2)
ax = symmetric_diff.plot(color="yellow")
closed_polygon_1.plot(ax=ax)
closed_polygon_2.plot(ax=ax)
