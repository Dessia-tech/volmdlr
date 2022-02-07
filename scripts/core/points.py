import volmdlr as vm

p1 = vm.Point3D(2.2409998778799998, 1.3700006035999999, -0.452188071297)
p2 = vm.Point3D(2.241, 1.37, -0.452187476119)
s = set([p1, p2])
if p1 == p2:
    print(p1.point_distance(p2))
    assert hash(p1) == hash(p2)