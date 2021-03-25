import math
import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw

for iarc, arc in enumerate([vme.Arc2D(vm.Point2D(0.3, 0.4),
                      vm.Point2D(0.35, 0.3),
                      vm.Point2D(0.37, 0.22)),
            vme.Arc2D(vm.Point2D(0.3, 0.4),
                      vm.Point2D(-0.35, 0.3),
                      vm.Point2D(0.37, 0.22))]):

    arc_triangle = vmw.ClosedPolygon2D([arc.start,arc.end, arc.center])

    ax = arc.plot()
    arc.center.plot(ax=ax)
    arc_triangle.plot(ax=ax)
    ax.set_aspect('equal')



    # Areas
    arc_area = arc.area()
    triangle_area = arc_triangle.area()
    arc_sl_area = arc.straight_line_area()
    if iarc == 0:
        area_diff = arc_area - triangle_area - abs(arc_sl_area)
    else:
        area_diff = arc_area + triangle_area - abs(arc_sl_area)

    assert math.isclose(area_diff, 0., abs_tol=1e-9)

    # COG
    arc_cog = arc.center_of_mass()
    triangle_cog = arc_triangle.center_of_mass()
    arc_sl_cog = arc.straight_line_center_of_mass()

    arc_cog.plot(ax=ax, color='r')
    triangle_cog.plot(ax=ax, color='g')
    arc_sl_cog.plot(ax=ax, color='b')

    if iarc == 0:
        cog_diff = (triangle_cog * triangle_area + arc_sl_cog * abs(
            arc_sl_area)) / (abs(arc_sl_area) + triangle_area) - arc_cog
    else:
        cog_diff = (triangle_cog * triangle_area + arc_cog * arc_area) / (arc_area + triangle_area) - arc_sl_cog
    print(cog_diff)

    assert math.isclose(cog_diff.x, 0., abs_tol=1e-9)
    assert math.isclose(cog_diff.y, 0., abs_tol=1e-9)


    # Second moment area
    arc_sma = arc.second_moment_area(arc.center)
    triangle_sma = arc_triangle.second_moment_area(arc.center)

    arc_slsma = arc.straight_line_second_moment_area(arc.center)

    if iarc == 0:
        Ix_diff = arc_sma[0] - triangle_sma[0] - abs(arc_slsma[0])
        Iy_diff = arc_sma[1] - triangle_sma[1] - abs(arc_slsma[1])
        Ixy_diff = arc_sma[2] - triangle_sma[2] - abs(arc_slsma[2])
    else:
        Ix_diff = arc_sma[0] + triangle_sma[0] - abs(arc_slsma[0])
        Iy_diff = arc_sma[1] + triangle_sma[1] - abs(arc_slsma[1])
        Ixy_diff = arc_sma[2] + triangle_sma[2] - abs(arc_slsma[2])

    print('triangle_sma', triangle_sma)
    print(Ix_diff, Iy_diff, Ixy_diff)
    assert math.isclose(Ix_diff, 0., abs_tol=1e-9)
    assert math.isclose(Iy_diff, 0., abs_tol=1e-9)
    # TODO: check this for arc > math.pi
    # assert math.isclose(Ixy_diff, 0., abs_tol=1e-9)

    # print(sma_arc, sma_triangle, sma_sl_arc)

