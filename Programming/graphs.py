import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import math

df = pd.read_csv('data.csv')

# Extract data from columns
time_values = df['Time']
velocity_values = df['Velocity']
acceleration_values = df['Acceleration']  # Assuming 'velocity' is the column for initial velocity

# Константы
g = 9.81  # ускорение свободного падения, м/с^2
v0 = 0.37  # начальная скорость, м/с
B = 1.29e-4  # коэффициент сопротивления воздуха
c = 0.045  # коэффициент формы
S = math.pi * 13 ** 2 /4 # площадь сечения, м^2
p0 = 1.225  # плотность воздуха на уровне моря, кг/м^3

# Другие параметры
E = (7.77e6*5 -6.77e6*5)/165# параметр, связанный с тягой
T = 1  # характерное время в экспоненте
m0 = 2900000  # начальная масса
a = 1300  # коэффициент в уравнении массы


# Система дифференциальных уравнений


def runge_kutta_4(f, t0, y0, x0, h0, t, *args):
    h = t - t0
    k1 = h * f(t0, y0, x0, h0, *args)[0]
    k2 = h * f(t0 + 0.5 * h, y0 + 0.5 * k1, x0, h0, *args)[0]
    k3 = h * f(t0 + 0.5 * h, y0 + 0.5 * k2, x0, h0, *args)[0]
    k4 = h * f(t0 + h, y0 + k3, x0, h0, *args)[0]
    k11 = h * f(t0, y0, x0, h0, *args)[1]
    k22 = h * f(t0 + 0.5 * h, y0, x0 + 0.5 * k11, h0, *args)[1]
    k33 = h * f(t0 + 0.5 * h, y0, x0 + 0.5 * k22, h0, *args)[1]
    k44 = h * f(t0 + h, y0, x0 + k33, h0, *args)[1]
    return y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6, x0 + (k11 + 2 * k22 + 2 * k33 + k44) / 6, y0


def solve_system(t, y, x, h, *args):
    vy = y
    vx = x
    h = h
    v= np.sqrt(vx**2 + vy**2)

    m = m0 - a * t
    Fthrust = 33832942*.7 + E * t
    a_t = np.degrees(90) * np.exp(-t / T)
    #print(Fthrust * np.sin(a_t) + Fthrust * np.cos(a_t), m * g, (1 / 2) * c * S * p0 * np.exp(-B * h) * v ** 2 * np.sin(a_t))
    # Your system of ODEs here

    return [
        (Fthrust * np.sin(a_t) - m * g  - (1 / 2) * c * S * p0 * np.exp(-B * h) * v ** 2 * np.sin(a_t)) / m,
        ((Fthrust * np.cos(a_t) - (1 / 2) * c * S * p0 * np.exp(-B * h) * vx ** 2 * np.cos(a_t)) / m) + v0 / t,
        vy
    ]


t0 = 1
y0 = [0, 0, 0]  # Initial conditions
t = np.linspace(t0, 100, 10000)  # Time grid
y = np.zeros((len(t), len(y0)))  # Solution grid
y[0] = y0

for i in range(1, len(t)):
    y[i] = runge_kutta_4(solve_system, t[i - 1], y[i - 1][0], y[i - 1][1], y[i - 1][2], t[i])

vx = [y[i][1] for i in range(len(y))]
vy = [y[i][2] for i in range(len(y))]
#print(S)

# Вычисляем сумму квадратов vx и vy и корень из этой суммы
speed_magnitude = [np.sqrt(vx[i] ** 2 + vy[i] ** 2) for i in range(len(vx))]

velocity_values_interp = np.interp(t, time_values, velocity_values)

# Построение графика
plt.figure(figsize=(10, 6))
plt.plot(t, speed_magnitude, label='Мат Модель')
plt.plot(t, velocity_values_interp, label='KSP')
plt.title('Скорость (модуль) от времени')
plt.xlabel('Время, с')
plt.ylabel('Скорость, м/с')
plt.legend()
plt.grid(True)
plt.show()