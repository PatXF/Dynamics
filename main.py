from pwe import *
from scd import *
from egi import *
from nbm import *
from enbm import *
from nbmnl import *

# information needed
wn = 6.283        # omega natural in rads / sec
m = 0.2533        # discrete mass of the single degree of freedom system in kip - sec^2 / in
k = 10            # stiffness in kips / in
zeta = 0.05       # damping ratio
P0 = 10           # loading amplitude
del_t = 0.1       # time steps in secs
t = 0.6           # time of application of load in secs
tol = "sin"       # type of load
toa = 1           # total time of response under study in secs
wl = np.pi / 0.6  # omega for load in rads / sec
ui = 0            # initial displacement
ui_ = 0           # initial velocity
gamma = 0.5       # gamma value for newark-beta method analysis
beta = 0.25       # beta value for newark-beta method analysis


def get_f(u, u_prev):  # definition of the stiffness function for non-linear analysis
    if u > u_prev:
        f_prev = 10 * u_prev
        if f_prev > 7.5:
            f_prev = 7.5
        f = f_prev + 10 * (u - u_prev)
        if f > 7.5:
            return 7.5
        else:
            return f
    else:
        f_prev = 10 * u_prev
        if f_prev < -7.5:
            f_prev = -7.5
        f = f_prev - 10 * (u_prev - u)
        if f < -7.5:
            return -7.5
        else:
            return f


solve_pwe(wn, m, k, zeta, P0, del_t, t, tol, toa, wl, ui, ui_)
solve_scd(wn, m, k, zeta, P0, del_t, t, tol, toa, wl, ui, ui_)
solve_egi(wn, m, k, zeta, P0, del_t, t, tol, toa, wl, ui, ui_)
solve_nbm(wn, m, k, zeta, P0, del_t, t, tol, toa, wl, ui, ui_, gamma, beta)
solve_enbm(wn, m, k, zeta, P0, del_t, t, tol, toa, wl, ui, ui_, gamma, beta)
solve_nbmnl(wn, m, k, zeta, P0, del_t, t, tol, toa, wl, ui, ui_, gamma, beta, get_f)