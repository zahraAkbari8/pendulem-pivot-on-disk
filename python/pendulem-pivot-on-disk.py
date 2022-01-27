#!/usr/bin/env python

import numpy as np
from numpy import sin, cos, pi, sqrt, floor
import matplotlib.pyplot as plt
from  scipy.integrate import solve_ivp

START_t = 0
STOP_t = 10*pi
EPSILON = 0.001
NUMBER = int((STOP_t - START_t) / EPSILON)
Y0 = [-pi/2, 0, pi/3, 0]

M = 1
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


t = np.linspace(START_t, STOP_t, NUMBER)
sol = solve_ivp(prime_func, [START_t, STOP_t], Y0, t_eval=t, method = "BDF")
phi = sol.y[0]
theta = sol.y[2]
x = a*cos(phi) + b*sin(theta)
y = a*sin(phi) - b*cos(theta)

fig = plt.figure()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(223)
ax3 = fig.add_subplot(122)
ax1.plot(sol.t, phi)
ax2.plot(sol.t, theta)
ax3.plot(x, y)

ax1.set_xlabel("$t$")
ax1.set_ylabel("$\\phi$")
ax2.set_xlabel("$t$")
ax2.set_ylabel("$\\theta$")
ax3.set_xlim([-(a+b+1), (a+b+1)])
ax3.set_ylim([-(a+b+1), (a+b+1)])
ax3.grid()

fig.set_size_inches((20, 11))
plt.savefig("fig: "+"[M, m, a, b] ="+str([M, m, a, b])+"-"+"Y0 ="+str(Y0)
            +"time = "+str(STOP_t - START_t)+".png")
plt.show()
