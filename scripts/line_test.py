import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
from volmdlr import plot_data


p1 = vm.Point2D([0,0])
p2 = vm.Point2D([1,1])
line = vm.Line2D(p1,p2)
contour = vm.Contour2D([line])
plot_datas = contour.plot_data()
sol = [plot_datas.to_dict()]
plot_data.plot_d3(sol)