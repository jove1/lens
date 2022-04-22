#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot  as plt

def cast(ray, lens):
    path = []
    x, y, dy = ray
    path.append( ray )
    for xl, sl in lens:
        y += xl*dy
        dy -= y*sl
        x += xl
        path.append( (x, y, dy) )
    return np.asarray(path)

from matplotlib.widgets import Slider, RadioButtons
from matplotlib.collections import LineCollection
x_sample = 0
x_screen = 8
x_l1 = 2
x_l2 = 6

fig = plt.figure()
w_mag = Slider(fig.add_axes([0.1,0.90,0.8,0.05]), 'mag', 1, 10)

def update(val=None):
    mag = w_mag.val

    x_img = x_l2 - (x_screen-x_l2)/mag
    s_l1 = 1/(x_img-x_l1) + 1/(x_l1-x_sample)
    s_l2 = 1/(x_l2-x_img) + 1/(x_screen-x_l2)
    print(x_img, s_l1, s_l2)

    paths = []
    for y in np.linspace(-1, 1, 2):
        for dy in [0, 
                -y/(x_l1-x_sample), 
                +y/(x_l1-x_sample), 
                #-y/(x_l1-x_sample-1/s_l1)
                ]:
#        for dy in np.linspace(-1, 1, 3):
            path = cast( (0, y, dy), [ (x_l1-x_sample, s_l1), (x_l2-x_l1, s_l2), (x_screen-x_l2, 0) ] )
            paths.append(path[:,:2])

    l.set_segments(paths)
    fig.canvas.draw_idle()

w_mag.on_changed(update)

fig = plt.figure()
ax = fig.gca()
l = LineCollection([], linewidths=0.5, colors='k')
ax.add_collection(l)
#ax.grid(0)
ax.axhline(0, ls=":", c="k")

# lens
ax.axvline(x_l2)
ax.axvline(x_l1)

# planes
ax.axvline(x_sample, ls="--")
ax.axvline(x_screen, ls="--")

ax.set_xlim(x_sample, x_screen)
ax.set_ylim(-20, 20)
update()

plt.show()
