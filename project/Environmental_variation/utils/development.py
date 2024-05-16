import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

time = np.linspace(0,100,100)

A = 0.9 
B = 0
L = 10 
R = 8
params = A, B, L, R
y0 = 0.0

def dENV_dt(y,t, params):

    A, B, L, R  = params
    epsilon = np.random.normal(0, 1)
    sin =  A * np.sin(2 * np.pi * t / (L * R)) + B * epsilon 
    print(sin)
    return sin

variation = odeint(dENV_dt, y0, time, args=(params,))

env = [dENV_dt(y, time, params) for time, y in zip(time, variation)]

plt.plot(time, env, linestyle="dashdot", color="blue", label="sin")
# plt.plot(time, variation, label="env inside odeint")
# # plt.plot(time,variation[:,1], linestyle="dashdot", color="red", label="cos")
# # plt.plot(time,variation[:,2], linestyle="dashdot", color="purple", label="sum")
plt.title("Environmental Variation")
plt.xlabel("Time")
plt.ylabel("Variation")
plt.grid()
plt.legend()
plt.show()
