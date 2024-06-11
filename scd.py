# Second central difference formulation to solve structural dynamics problem

import numpy as np
import matplotlib.pyplot as plt


def oscillate(ot, theta):
    if ot == "sin":
        return np.sin(theta)
    elif ot == "cos":
        return np.cos(theta)


def solve_scd(wn, m, k, zeta, P0, del_t, t, tol, toa, wl, ui, ui_):

    # initial calculations
    steps = int(toa / del_t)
    c = zeta * 2 * m * wn
    p_ini = P0 * oscillate(tol, 0)
    ui__ = (p_ini - c * ui_ + k * ui) / m
    uim1 = ui - del_t * ui_ + ((del_t ** 2) / 2) * ui__
    disps = [ui]
    velo = [ui_]
    times = [0]
    k_hat = m / (del_t ** 2) + c / (2 * del_t)

    for i in range(steps):
        toc = del_t * (i + 1)  # time of consideration for the current step
        pi = 0
        if toc <= t:
            pi = np.round(P0 * oscillate(tol, wl * toc), 4)
        p_hat = pi - ((m / (del_t ** 2)) - (c / (2 * del_t))) * uim1 - (k - (2 * m) / (del_t ** 2)) * ui
        ui1 = np.round(p_hat / k_hat, 4)
        ui1_ = np.round((ui1 - uim1) / (2 * del_t), 4)
        disps.append(ui1)
        velo.append(ui1_)
        uim1 = ui
        ui = ui1
        ui_ = ui1_
        times.append(toc)

    data_obtained_in_book = [0, 0.1914, 0.6293, 1.1825, 1.5808, 1.5412, 0.9141, -0.0247, -0.8968, -1.3726, -1.2940]
    theoretical_disps = [0, 0.0328, 0.2332, 0.6487, 1.1605, 1.5241, 1.4814, 0.9245, 0.0593, -0.7751, -1.2718]
    plt.figure(figsize=(10, 6))
    plt.plot(times, disps, label='Displacement (x)')
    plt.scatter(times, theoretical_disps, color="red", label="theoretical displacement")
    plt.scatter(times, data_obtained_in_book, color="green", label="data obtained in book")
    plt.title('SDOF System Response using Second Central Difference Method vs theoretical results vs data obtained in book')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [in]')
    plt.legend()
    plt.grid()
    plt.show()