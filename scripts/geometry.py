import volmdlr.geometry as vmg
import math
import random

for (x1, x2, cos_min, cos_max) in [(0., 0.1, math.cos(0.1), 1),
                                   (0.1, 0.2, math.cos(0.2), math.cos(0.1))]:
    fx1, fx2 = vmg.cos_image(x1, x2)
    assert fx1 == cos_min
    assert fx2 == cos_max


# Visual test
x = [i/200.*4*math.pi-2*math.pi for i in range(200)]
cosx = [math.cos(xi) for xi in x]
sinx = [math.sin(xi) for xi in x]

import matplotlib.pyplot as plt



for i in range(3):
    fig, (ax1, ax2) = plt.subplots(2)
    ax1.plot(x, cosx)
    ax2.plot(x, sinx)

    x1 = random.uniform(-2*math.pi, 2*math.pi)
    x2 = random.uniform(-2 * math.pi, 2 * math.pi)
    x1, x2 = sorted([x1,x2])
    cos_min, cos_max = vmg.cos_image(x1, x2)
    sin_min, sin_max = vmg.sin_image(x1, x2)

    ax1.plot([x1, x2, x2, x1, x1], [cos_min, cos_min, cos_max, cos_max, cos_min], color='r')
    ax2.plot([x1, x2, x2, x1, x1], [sin_min, sin_min, sin_max, sin_max, sin_min], color='r')
    # print(x1, x2)
