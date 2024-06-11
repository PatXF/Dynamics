##code to analyse a structural dynamics problem using piecewise exact method (single degree of freedom)
import numpy as np
import matplotlib.pyplot as plt


def oscillate(ot, theta):
    if ot == "sin":
        return np.sin(theta)
    elif ot == "cos":
        return np.cos(theta)


def solve_pwe(wn, m, k, zeta, P0, del_t, t, tol, toa, wl, ui, ui_):

    # initial calculations
    E = np.round(np.exp(-(zeta * wn * del_t)), 4)  # checked with example, working fine
    wd = np.round(wn * np.sqrt(1 - (zeta ** 2)), 4)  # checked with example, working fine
    swd = np.round(np.sin(wd * del_t), 4)  # checked with example, working fine
    cwd = np.round(np.cos(wd * del_t), 4)  # checked with example, working fine

    A = np.round(E * ((zeta / np.sqrt(1 - (zeta ** 2))) * swd + cwd), 4)  # checked with example, working fine
    B = np.round(E * (1 / wd) * swd, 4)  # checked with example, working fine
    C = np.round((((2 * zeta) / (wn * del_t)) + E * (((1 - 2 * (zeta ** 2)) / (wd * del_t) - (zeta / np.sqrt(1 - (zeta ** 2)))) * swd - (1 + (2 * zeta) / (wn * del_t)) * cwd)) / k, 4)  # checked with example, working fine
    D = np.round((1 - (2 * zeta) / (wn * del_t) + E * (((2 * (zeta ** 2) - 1) / (wd * del_t)) * swd + ((2 * zeta) / (wn * del_t)) * cwd)) / k, 4)  # checked with example, working fine
    A_ = np.round(-E * ((wn * swd) / np.sqrt(1 - (zeta ** 2))), 4)  # checked with example, working fine
    B_ = np.round(E * (cwd - (zeta * swd) / np.round(1 - (zeta ** 2))), 4)  # checked with example, working fine
    C_ = np.round((((-1 / del_t) + E * (swd * ((wn / np.sqrt(1 - zeta ** 2)) + zeta / (del_t * np.sqrt(1 - zeta ** 2))) + cwd / del_t)) / k), 4)  # checked with example, working fine
    D_ = np.round((1 / (k * del_t)) * (1 - E * ((zeta * swd) / np.sqrt(1 - zeta ** 2) + cwd)), 4)  # checked with example, working fine

    total_steps = int(toa / del_t)

    disps = [ui]
    velo = [ui_]
    times = [0]

    for i in range(0, total_steps):
        toc = del_t * i  # time of consideration for the current step
        toc1 = del_t * (i + 1)  # time of consideration for the next step
        pi = 0
        pi1 = 0
        if toc <= t:
            pi = np.round(P0 * oscillate(tol, wl * toc), 4)
        if toc1 <= t:
            pi1 = np.round(P0 * oscillate(tol, wl * toc1), 4)
        ui1 = np.round(A * ui + B * ui_ + C * pi + D * pi1, 4)
        ui1_ = np.round(A_ * ui + B_ * ui_ + C_ * pi + D_ * pi1, 4)
        disps.append(ui1)
        velo.append(ui1_)
        ui = ui1
        ui_ = ui1_
        times.append(toc1)

    theoretical_disps = [0, 0.0328, 0.2332, 0.6487, 1.1605, 1.5241, 1.4814, 0.9245, 0.0593, -0.7751, -1.2718]
    plt.figure(figsize=(10, 6))
    plt.plot(times, disps, label='Displacement (x)')
    plt.scatter(times, theoretical_disps, color="red", label="theoretical displacements")
    plt.title('SDOF System Response using Piecewise Exact Method vs theoretical results')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [in]')
    plt.legend()
    plt.grid()
    plt.show()