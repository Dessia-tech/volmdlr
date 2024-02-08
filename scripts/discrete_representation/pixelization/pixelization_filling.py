"""
Showcase of filling a pixelization.
"""

from volmdlr import Point2D
from volmdlr.core import EdgeStyle
from volmdlr.discrete import MatrixBasedPixelization
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

closed_polygon = ClosedPolygon2D(points)
pixelization = MatrixBasedPixelization.from_closed_polygon(closed_polygon, PIXEL_SIZE)

# Filling the outer pixels
outer_filled_pixelization = pixelization.fill_outer_pixels()
ax = outer_filled_pixelization.plot()
closed_polygon.plot(ax=ax, edge_style=EdgeStyle(color="r"))

# Filliong the inner pixels
inner_filled_pixelization = pixelization.fill_enclosed_pixels()
ax = inner_filled_pixelization.plot()
closed_polygon.plot(ax=ax, edge_style=EdgeStyle(color="r"))
