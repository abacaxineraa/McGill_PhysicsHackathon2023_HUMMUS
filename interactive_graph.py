import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import scipy.ndimage.filters as filters

t = 0
n = 1000
l = 10
c = 1
damping = 0.01
h = 5
d = 1

def update(frame):
    global t
    t += 0.5
    y = np.zeros_like(x)
    for j in range(1, n):
        An = (2 * h * l ** 2) / ((j ** 2) * (math.pi ** 2) * d * (l * d)) * np.sin(d * j * math.pi / l)
        y += An * np.sin(j * math.pi * x / l) * np.cos(2 * np.pi * j * c * t / (2 * l)) * math.e ** (-damping * t)
    y = filters.gaussian_filter1d(y, sigma=30)
    line.set_ydata(y)
    return line,

x = np.linspace(0, l, 1000)
y = np.zeros_like(x)
fig, ax = plt.subplots()
line, = ax.plot(x, y)
ax.set_xlim([0, l])
ax.set_ylim(-10, 10)

def onclick(event):
    global h, d, ani, t
    h = event.ydata*2
    d = event.xdata
    ani.event_source.stop()
    t = 0
    ani = FuncAnimation(fig, update, frames=200, blit=True, interval=10)
    plt.show()
    print(h)
    print(d)

cid = fig.canvas.mpl_connect('button_press_event', onclick)

ani = FuncAnimation(fig, update, frames=200, blit=True, interval=10)
plt.show()
