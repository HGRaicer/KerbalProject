import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import math

df = pd.read_csv('data.csv')

# Извлечение данных из столбцов
time_values = df['Time']
velocity_values = df['Velocity']
acceleration_values = df['Acceleration']

# Константы
g = 9.81  # ускорение свободного падения, м/с^2
v0 = 0.37  # начальная скорость, м/с
B = 1.29e-4  # коэффициент сопротивления воздуха
c = 0.045  # коэффициент лобового сопротивления
S = math.pi * 13 ** 2 / 4  # площадь сечения, м^2
p0 = 1.225  # плотность воздуха на уровне моря, кг/м^3
# Другие параметры
E = (7.77e6 * 5 - 6.77e6 * 5) / 165  # коэффициент изменение тяги
T = 170  # характерное время в экспоненте
m0 = 2900000  # начальная масса
a = 1300  # коэффициент сжигания топлива
KSP = 1/320 # Коэффициент сопостовления КСП и мат модели


# Система дифференциальных уравнений


def runge_kutta_4(f, f2, f3, t0, vx0, vy0, h0, t, *args):
    h = t - t0
    k11 = f(t0, vx0, vy0, h0, *args)
    k12 = f2(t0, vx0, vy0, h0, *args)
    k13 = f3(t0, vx0, vy0, h0, *args)
    k21 = f(t0 + 0.5 * h, vx0 + 0.5 * k11 * h, vy0 + 0.5 * k12 * h, h0 + k13 * h / 2, *args)
    k22 = f2(t0 + 0.5 * h, vx0 + 0.5 * k11 * h, vy0 + 0.5 * k12 * h, h0 + k13 * h / 2, *args)
    k23 = f3(t0 + 0.5 * h, vx0 + 0.5 * k11 * h, vy0 + 0.5 * k12 * h, h0 + k13 * h / 2, *args)
    k31 = f(t0 + 0.5 * h, vx0 + 0.5 * k21 * h, vy0 + 0.5 * k22 * h, h0 + k23 * h / 2, *args)
    k32 = f2(t0 + 0.5 * h, vx0 + 0.5 * k21 * h, vy0 + 0.5 * k22 * h, h0 + k23 * h / 2, *args)
    k33 = f3(t0 + 0.5 * h, vx0 + 0.5 * k21 * h, vy0 + 0.5 * k22 * h, h0 + k23 * h / 2, *args)
    k41 = f(t0 + 0.5 * h, vx0 + k31 * h, vy0 + 0.5 * k32 * h, h0 + k33 * h, *args)
    k42 = f2(t0 + 0.5 * h, vx0 + k31 * h, vy0 + 0.5 * k32 * h, h0 + k33 * h, *args)
    k43 = f3(t0 + 0.5 * h, vx0 + k31 * h, vy0 + 0.5 * k32 * h, h0 + k33 * h, *args)
    return vx0 + (k11 + 2 * k21 + 2 * k31 + k41) * h / 6 + v0 / t, vy0 + (k12 + 2 * k22 + 2 * k32 + k42) * h / 6, h0+(
                k13 + 2 * k23 + 2 * k33 + k43) * h / 6


def solve_system_vy(t, vx, vy, h, *args):
    v = np.sqrt(vx ** 2 + vy ** 2)
    m = m0 - a * t
    Fthrust = 6.77e6 * 5 + E * t
    a_t = np.degrees(90) * np.exp(-t / T)

    return (abs(Fthrust * np.sin(a_t)) - m * g - (1 / 2) * c * S * p0 * np.exp(-B * h) * v ** 2 * np.sin(a_t)) / m


def solve_system_vx(t, vx, vy, h, *args):

    v = np.sqrt(vx ** 2 + vy ** 2)
    m = m0 - a * t
    Fthrust = 6.77e6 * 5 + E * t
    a_t = np.degrees(90) * np.exp(-t / T)

    return ((abs(Fthrust * np.cos(a_t)) - (1 / 2) * c * S * p0 * p0 * np.exp(-B * h) * v ** 2 * np.cos(a_t)) / m) +v0/t


def change_h(t, vx, vy,h, *args):
    return abs(vy)


t0 = 1  # Первоначальные условия
num = 1200
stop = 120
t = np.linspace(t0, stop, num)  # Временная сетка
ans = [[0,0,0]]
sec = stop / num
for i in range(1, num):
    time = i*sec
    ans.append(
        runge_kutta_4(solve_system_vx, solve_system_vy, change_h, t0, ans[i - 1][0], ans[i - 1][1], ans[i - 1][2],
                      time))

hh = []
speed_magnitude = []
for idx in range(len(ans)):
    speed_magnitude.append((math.sqrt(ans[idx][0] ** 2 + ans[idx][1] ** 2))*KSP)
    hh.append(ans[idx][2])

velocity_values_interp = np.interp(t, time_values, velocity_values)

tt = np.linspace(t0, len(speed_magnitude), len(speed_magnitude))

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
