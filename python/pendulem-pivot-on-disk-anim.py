#!/usr/bin/env python

import numpy as np
from numpy import sin, cos, pi, sqrt, floor
from  time import time
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.patches as mpatches
from  scipy.integrate import solve_ivp


START_t = 0
STOP_t = 30*pi
EPSILON = 0.001
NUMBER = int((STOP_t - START_t) / EPSILON)
ANIM_SPEED = 250
FRAMES = int(NUMBER/ANIM_SPEED)
# FRAMES = 100
Y0 = [-pi/2+0.1, 0, pi+0.1, 0]

M = sqrt(1717)
m = 1
a = 1
b = sqrt(17)
g = 9.78

def prime_func(t, Y):
    phi = Y[0]
    p_phi = Y[1]
    theta = Y[2]
    p_theta = Y[3]
    pp_phi = m * cos ( theta - phi ) / ( M + m * cos ( theta -phi ) **2 ) * ( b * ( ( p_theta )**2 ) - a * ( ( p_phi ) **2 ) * sin(theta - phi))
    pp_phi += m*g*(sin(theta) *sin(theta - phi) - cos(phi)) / (a*(M+m*cos(theta -phi)**2))
    pp_theta = (cos(theta - phi)) / (M+m*(cos(theta-phi)**2)) * ((M+m)*a*(p_phi**2)-m*(p_theta**2)*sin(theta-phi))
    pp_theta += (g/(b*(M+m*cos(theta-phi)**2)))*(m*(cos(phi))*(sin(theta-phi)) - (M+m)*sin(theta))

    return [p_phi, pp_phi, p_theta, pp_theta]


def anim_init():
    ax.set_xlim([-(a+b+1), (a+b+1)])
    ax.set_ylim([-(a+b+1), (a+b+1)])
    ax.set_aspect("equal")
    ax.grid()
    return line, line2, circle, end

def anim_update(frame):
    line.set_data(x[0:frame*ANIM_SPEED], y[0:frame*ANIM_SPEED])
    line2.set_data([xc[frame*ANIM_SPEED], x[frame*ANIM_SPEED]], [yc[frame*ANIM_SPEED], y[frame*ANIM_SPEED]])
    circle.set_data(xc[frame*ANIM_SPEED], yc[frame*ANIM_SPEED])
    end.set_data(x[frame*ANIM_SPEED], y[frame*ANIM_SPEED])
    return line, line2, circle, end

t = np.linspace(START_t, STOP_t, NUMBER)
sol = solve_ivp(prime_func, [START_t, STOP_t], Y0, t_eval=t, method = "BDF")
print(sol.message)

phi = sol.y[0]
theta = sol.y[2]
x = a*cos(phi) + b*sin(theta)
y = a*sin(phi) - b*cos(theta)
xc = a*cos(phi)
yc = a*sin(phi)

fig = plt.figure()
ax = fig.add_subplot(111)

line, = ax.plot([], [], color="blue")
line2, = ax.plot([],[], color="black")
circle = ax.add_patch(mpatches.Circle((0,0), radius = a, fill = False))
circle, = ax.plot([],[], "bo")
end, = ax.plot([],[], "ro")

starting_time = time()
ani = FuncAnimation(fig, func=anim_update, frames=np.arange(0,FRAMES),
                    init_func = anim_init, blit = True, repeat =False)
animation_time = time()
print("animation done in ." + str(animation_time - starting_time))
ani.save("anim.gif", writer="ffmpeg", fps=24, dpi=300)
print("animation saved in ."+str(time() - animation_time))
# fig.set_size_inches((20, 11))
# plt.show()
print("Done!")
