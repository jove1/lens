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

x_cl = 5
x_sample = 10
x_obj = 14
x_diff = 17
x_proj = 30
x_screen = 40

s_ol = 1/(x_diff-x_obj)
x_img = x_obj + 1/(s_ol - 1/(x_obj-x_sample))

from matplotlib.widgets import Slider, RadioButtons
from matplotlib.collections import LineCollection

fig = plt.figure()
y = 1
y -= 0.1; w_shift = Slider(fig.add_axes([0.1,y,0.8,0.05]), 'Shift', -20, 20, 0)
y -= 0.1; w_cl = Slider(fig.add_axes([0.1,y,0.8,0.05]), 'CL', 0, 1, 0)
y -= 0.2; w_proj = RadioButtons(fig.add_axes([0.1,y,0.12,0.15]), ['IMG', 'DIFF']); 

def update(val=None):
    s_cl = 1/x_cl + w_cl.val/(x_sample-x_cl)
    if w_proj.value_selected == 'IMG':
        s_proj = 1/(x_proj-x_img) + 1/(x_screen-x_proj)
    else: # 'DIFF'
        s_proj = 1/(x_proj-x_diff) + 1/(x_screen-x_proj)

    paths = []
    for dy in np.linspace(-1, 1, 2):
        path = cast( (0, w_shift.val, dy), [ (x_cl, s_cl), (x_sample-x_cl, 0) ] ) # condensor
        paths.append(path[:,:2])
        for dy in np.linspace(-4, 4, 3): # diffraction on sample
            x0, y0, dy0 = path[-1]
            path2 = cast( (x0,y0,dy0+dy), [ (x_obj-x_sample, s_ol), (x_proj-x_obj, s_proj), (x_screen-x_proj, 0) ] )
            paths.append(path2[:,:2])

    l.set_segments(paths)
    fig.canvas.draw_idle()

w_cl.on_changed(update)
w_shift.on_changed(update)
w_proj.on_clicked(update)

fig = plt.figure()
ax = fig.gca()
l = LineCollection([], linewidths=0.5, colors='k')
ax.add_collection(l)
ax.grid(0)
ax.axhline(0, ls=":", c="k")

# lens
ax.axvline(x_cl)
ax.axvline(x_obj)
ax.axvline(x_proj)

# planes
ax.axvline(x_sample, ls="--")
ax.axvline(x_screen, ls="--")
ax.axvline(x_diff, ls="--")
ax.axvline(x_img, ls="--")

ax.set_xlim(0, x_screen)
ax.set_ylim(-50,50)
update()

plt.show()
