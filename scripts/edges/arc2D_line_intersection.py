import volmdlr
import volmdlr.edges
import matplotlib.pyplot as plt

# Random arc
# i = volmdlr.Point2D.random(-1,1,-1,1)
# e = volmdlr.Point2D.random(-1,1,-1,1)
# s = volmdlr.Point2D.random(-1,1,-1,1)



i = volmdlr.Point2D(-0.2, -0.25)
e = volmdlr.Point2D(-1, 2)
s = volmdlr.Point2D(1.5, 2)

print('i', i)
print('s', s)
print('e', e)

a = volmdlr.edges.Arc2D(s, i, e)
start1 = volmdlr.Point2D(-2, -1.5)
end1 = volmdlr.Point2D(1.25, 1.75)
start2 = volmdlr.Point2D(1, 0)
end2 = volmdlr.Point2D(2, -2)
start3 = volmdlr.Point2D(0.7, 1.8)
end3 = volmdlr.Point2D(2, 2.2)
l1 = volmdlr.edges.LineSegment2D(start1, end1)
l2 = volmdlr.edges.LineSegment2D(start2, end2)
l3 = volmdlr.edges.LineSegment2D(start3, end3)


ax = a.plot()
l1.plot(ax=ax, color='r')
l2.plot(ax=ax, color='b')
l3.plot(ax=ax, color='g')
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)


point1 = a.linesegment_intersections(l1)
point2 = a.linesegment_intersections(l2)
point3 = a.linesegment_intersections(l3)
# print('l1', point1)
# print('l2', point2)
# print('l3', point3)

assert len(point1) == 1
assert len(point2) == 1
assert len(point3) == 0
# l1 = volmdlr.edges.LineSegment2D

i = volmdlr.Point2D(-0.2, -0.25)
e = volmdlr.Point2D(1, 2)
s = volmdlr.Point2D(-1.5, -2)

a = volmdlr.edges.Arc2D(s, i, e)
start1 = volmdlr.Point2D(-2, -1.5)
end1 = volmdlr.Point2D(1.25, 1.75)
start2 = volmdlr.Point2D(1, 0)
end2 = volmdlr.Point2D(2, -2)
start3 = volmdlr.Point2D(0.7, 1.8)
end3 = volmdlr.Point2D(2, 2.2)
l1 = volmdlr.edges.LineSegment2D(start1, end1)
l2 = volmdlr.edges.LineSegment2D(start2, end2)
l3 = volmdlr.edges.LineSegment2D(start3, end3)


ax = a.plot()
l1.plot(ax=ax, color='r')
l2.plot(ax=ax, color='b')
l3.plot(ax=ax, color='g')
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)


point1 = a.linesegment_intersections(l1)
point2 = a.linesegment_intersections(l2)
point3 = a.linesegment_intersections(l3)
# print('l1', point1)
# print('l2', point2)
# print('l3', point3)

assert len(point1) == 1
assert len(point2) == 0
assert len(point3) == 1
