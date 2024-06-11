# Code to perform analysis on a structural dynamics problem using Newmark beta method (explicit)

import numpy as np
import matplotlib.pyplot as plt


def oscillate(ot, theta):
    if ot == "sin":
        return np.sin(theta)
    elif ot == "cos":
        return np.cos(theta)


def solve_enbm(wn, m, k, zeta, P0, del_t, t, tol, toa, wl, ui, ui_, y, b):

    # additional information obtained from information provided
    cc = 2 * m * wn
    c = zeta * cc
    p_ini = P0 * oscillate(tol, 0)
    ui__ = (p_ini - c * ui_ - k * ui) / m
    steps = int(toa / del_t)
    times = [0]
    disps = [ui]
    velo = [ui_]
    k_hat = np.round(k + (y * c) / (b * del_t) + m / (b * (del_t ** 2)), 4)
    A = np.round(m / (b * del_t) + (y * c) / b, 4)
    B = np.round(m / (2 * b) + (c * del_t) * (y / (2 * b) - 1), 4)

    for i in range(steps):
        toc = i * del_t
        toc1 = (i + 1) * del_t
        times.append(toc1)
        p0 = P0 * oscillate(tol, wl * toc)
        if toc > t:
            p0 = 0
        p1 = P0 * oscillate(tol, wl * toc1)
        if toc1 > t:
            p1 = 0

        del_p = p1 - p0
        del_p_hat = del_p + A * ui_ + B * ui__
        del_u = del_p_hat / k_hat
        del_u_ = (y / (b * del_t)) * del_u - (y * ui_) / b + ui__ * del_t * (1 - y / (2 * b))
        del_u__ = del_u / (b * (del_t ** 2)) - ui_ / (b * del_t) - ui__ / (2 * b)
        ui = np.round(ui + del_u, 4)
        ui_ = np.round(ui_ + del_u_, 4)
        ui__ = np.round(ui__ + del_u__, 4)
        disps.append(ui)
        velo.append(ui_)

    theoretical_disps = [0, 0.0328, 0.2332, 0.6487, 1.1605, 1.5241, 1.4814, 0.9245, 0.0593, -0.7751, -1.2718]
    data_obtained_in_book = [0, 0.0437, 0.2326, 0.6121, 1.0825, 1.4309, 1.4231, 0.9622, 0.1908, -0.6044, -1.1442]
    plt.figure(figsize=(10, 6))
    plt.plot(times, disps, label='Displacement (x)')
    plt.scatter(times, theoretical_disps, color="red", label="theoretical displacement")
    plt.scatter(times, data_obtained_in_book, color="green", label="data obtained in book")
    plt.title('SDOF System Response using Newark-Beta Method (Explicit) \n vs theoretical results vs data obtained in book')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [in]')
    plt.legend()
    plt.grid()
    plt.show()