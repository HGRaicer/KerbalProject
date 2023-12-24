import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp


g = 1.622
J = 2980
P = 46695.6

y0 = 7500
x0 = 0
V0 = 550
theta0 = np.deg2rad(90)

t_span = (1, 100)
t_eval = np.linspace(*t_span, 500)


def mass_function(t):
    return 16400 - (P / J) * t


def odes(t, y):
    y_pos, x_pos, velocity, angle = y
    mass = mass_function(t)

    dydt = [-velocity * np.cos(angle),  # y_dot
            velocity * np.sin(angle),  # x_dot
            -P / mass + g * np.cos(angle),  # V_dot
            -g * np.sin(angle) / velocity]  # theta_dot
    return dydt


initial_values = [y0, x0, V0, theta0]

t_eval = np.linspace(t_span[0], t_span[1], 1000)

solution = solve_ivp(odes, t_span, initial_values, t_eval=t_eval)

print(solution.y)
# print(solution.y[1])
# print(solution.y[2])
# print(solution.y[3])


plt.title("Луна")
plt.xlabel("x")
plt.ylabel("y")
plt.plot(solution.y[1], solution.y[0])
plt.grid(True)
# plt.plot(x, y1, color="blue")
plt.show()
