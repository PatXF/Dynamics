# Code to perform non-linear analysis on a structural dynamics problem using Newmark beta method (explicit)

import numpy as np
import matplotlib.pyplot as plt


def oscillate(ot, theta):
    if ot == "sin":
        return np.sin(theta)
    elif ot == "cos":
        return np.cos(theta)


def nr(ui, fsi, del_pi, ki, a, del_t, get_f):
    fs1 = fsi
    del_r = del_pi
    k_hat = ki
    sum_u = 0
    while 1:
        del_u = del_r / k_hat
        u_temp = ui + del_u
        fsi1 = get_f(u_temp, ui)
        del_f = fsi1 - fs1 + (a / del_t) * del_u
        del_r = del_r - del_f
        fs1 = fsi1
        ui = u_temp
        sum_u += del_u
        if del_u <= 0.001:
            return np.round(ui, 4)


def solve_nbmnl(wn, m, k, zeta, P0, del_t, t, tol, toa, wl, ui, ui_, y, b, get_f):

    # initial calculations
    steps = int(toa / del_t)
    c = zeta * 2 * m * wn
    p_ini = P0 * oscillate(tol, 0)
    ui__ = (p_ini - c * ui_ + k * ui) / m  # checked with example, working fine
    A = np.round(m / (b * del_t) + (y * c) / b, 4)  # checked with example, working fine
    B = np.round(m / (2 * b) + del_t * (y / (2 * b) - 1) * c, 4)  # checked with example, working fine
    disps = [ui]
    velo = [ui_]
    times = [0]
    ui_prev = 0

    for i in range(steps - 1):
        toc = del_t * i  # time of consideration for the current step
        toc1 = del_t * (i + 1)  # time of consideration for the next step
        pi = 0
        pi1 = 0
        if toc <= t:
            pi = np.round(P0 * oscillate(tol, wl * toc), 4)
        if toc1 <= t:
            pi1 = np.round(P0 * oscillate(tol, wl * toc1), 4)
        del_p_hat = pi1 - pi + A * ui_ + B * ui__
        k = 0
        if ui != ui_prev:
            k = 10
        k_hat = k + (y * c) / (b * del_t) + (m / (b * (del_t ** 2)))
        fsi = get_f(ui, ui_prev)
        ui1 = nr(ui, fsi, del_p_hat, k_hat, A, del_t, get_f)
        del_u = ui1 - ui
        del_u_ = (y * del_u) / (b * del_t) - (y * ui_) / b + del_t * ui__ * (1 - y / (2 * b))
        del_u__ = (1 / (b * (del_t ** 2))) * del_u - ui_ / (b * del_t) - ui__ / (2 * b)
        ui = ui1
        ui_ = np.round(ui_ + del_u_, 4)
        ui__ = np.round(ui__ + del_u__, 4)
        disps.append(ui)
        velo.append(ui_)
        times.append(toc1)

    theoretical_disps = [0, 0.0437, 0.2326, 0.6121, 1.1143, 1.6213, 1.9889, 2.0947, 1.9233, 1.5593]
    plt.figure(figsize=(10, 6))
    plt.plot(times, disps, label='Displacement (x)')
    plt.scatter(times, theoretical_disps, color="red", label="theoretical displacements")
    plt.title('SDOF System Response using Newark-Beta method with modified Newton Raphson vs theoretical results')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [in]')
    plt.legend()
    plt.grid()
    plt.show()
