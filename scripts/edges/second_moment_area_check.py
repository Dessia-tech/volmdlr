import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw

arc = vme.Arc2D(vm.Point2D(0.3, 0.4),
                vm.Point2D(0.35, 0.3),
                vm.Point2D(0.37, 0.22))

# ax = arc.plot()
# arc.center.plot(ax=ax)
# ax.set_aspect('equal')

sma_arc = arc.second_moment_area(arc.center)
arc_triangle = vmw.ClosedPolygon2D([arc.start,arc.end, arc.center])
sma_triangle = arc_triangle.second_moment_area(arc.center)

sma_sl_arc = arc.straight_line_second_moment_area(arc.center)
Ix_diff = sma_arc[0] -  sma_triangle[0] - abs(sma_sl_arc[0])
Iy_diff = sma_arc[1] -  sma_triangle[1] - abs(sma_sl_arc[1])
Ixy_diff = sma_arc[2] -  sma_triangle[2] - abs(sma_sl_arc[2])
print('sma_triangle', sma_triangle)
print(sma_arc[0], sma_triangle[0], sma_sl_arc[0])
print(Ix_diff, Iy_diff, Ixy_diff)
assert Ix_diff == 0.
assert Iy_diff == 0.
assert Ixy_diff == 0.

print(sma_arc, sma_triangle, sma_sl_arc)