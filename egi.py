# Code for solving a structural dynamics problem using euler gauss integral method

import numpy as np
import matplotlib.pyplot as plt


def oscillate(ot, theta):
    if ot == "sin":
        return np.sin(theta)
    elif ot == "cos":
        return np.cos(theta)


def calc_acceleration(m, p0, c, v_0, k, v0):
    return (1 / m) * (p0 - c * v_0 - k * v0)


def calc_velocity(v_0, t, acc0, acc1):
    return v_0 + (t / 2) * (acc0 + acc1)


def calc_displacement(v0, v_0, t, acc0, acc1):
    return v0 + v_0 * t + ((t ** 2) / 4) * (acc0 + acc1)


def solve_egi(wn, m, k, zeta, P0, del_t, t, tol, toa, wl, ui, ui_):

    # additional information obtained from information provided
    cc = 2 * m * wn
    c = zeta * cc
    steps = int(toa / del_t)
    disps = [ui]
    velo = [ui_]
    times = [0]

    for i in range(steps):
        toc = i * del_t
        toc1 = (i + 1) * del_t
        times.append(toc1)
        p0 = P0 * oscillate(tol, wl * toc)
        if toc > t:
            p0 = 0
        acc0 = calc_acceleration(m, p0, c, ui_, k, ui)
        acc1 = 0
        while 1:
            v_1 = calc_velocity(ui_, del_t, acc0, acc1)
            v1 = calc_displacement(ui, ui_, del_t, acc0, acc1)
            p1 = P0 * oscillate(tol, wl * toc1)
            if toc1 > t:
                p1 = 0
            temp_acc1 = calc_acceleration(m, p1, c, v_1, k, v1)
            if abs(temp_acc1 - acc1) <= 0.001:
                ui = np.round(v1, 4)
                ui_ = np.round(v_1, 4)
                disps.append(ui)
                velo.append(ui_)
                break
            acc1 = temp_acc1

    theoretical_disps = [0, 0.0328, 0.2332, 0.6487, 1.1605, 1.5241, 1.4814, 0.9245, 0.0593, -0.7751, -1.2718]
    plt.figure(figsize=(10, 6))
    plt.plot(times, disps, label='Displacement (x)')
    plt.scatter(times, theoretical_disps, color="red", label="theoretical_disps")
    plt.title('SDOF System Response using Euler Gauss Integral vs theoretical results')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [in]')
    plt.legend()
    plt.grid()
    plt.show()