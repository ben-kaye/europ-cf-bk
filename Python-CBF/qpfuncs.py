import numpy as np
import scipy as sp
from scipy import sparse

def sim_step(x, u, step_sz):
    phi = x[2]

    Ax = np.array(((np.cos(phi), 0),(np.sin(phi), 0),(0, 1)))

    return x + step_sz*Ax.dot(u)

def ref_step(r, step_sz):
    phi_r = r[2]

    r_dot = np.array(((np.cos(phi_r), 0), (np.sin(phi_r), 0), (0, 1)))
    return r + step_sz*np.hstack([r_dot.dot(r[[4,3]]), 0, r[5], 0])  
    
def saturate_ctrls(u, max_turn, v_min, v_max):
    v = u[0]
    v = max(v, v_min)
    v = min(v, v_max)

    omega = u[1]
    omega = np.sign(omega)*max(abs(omega), max_turn)

    return np.array((v, omega))

def get_sig_bf_constr(x, u, p_o, max_turn, delta):
    p = x[:2]
    phi = x[2]
    v = u[0]

    dir_v = np.array((np.cos(phi), np.sin(phi)))

    p_xo = p_o - p

    cross_k = np.array(((0, 1),(-1, 0))) 

    sign_r = np.sign(dir_v.dot(cross_k.dot(p_xo)))
    R_min = v/max_turn

    r = sign_r*R_min*cross_k.dot(dir_v)

    z = p_xo - r

    h = z.dot(z) - (delta + R_min)**2

    q = np.exp(1/h**2)

    sig = (q - 1)/(q + 1)

    # Abf = np.array((0, -1))
    ubf = sig*sign_r*max_turn

    return ubf

def clf_controls(x, r, k1, k2):
    p = x[:2]
    phi = x[2]
    p_r = r[:2]
    phi_r = r[2]
    phi_rdot = r[3]
    v_r = r[4]
    v_rdot = r[5]

    e3 = phi_r - phi
    c = np.cos(phi)
    s = np.sin(phi)
    e1, e2 = np.array(((c, s), (-s, c))).dot(p_r - p)

    alpha = np.arctan(e2/v_r)
    e3_aux = e3 + alpha

    fr = (v_r**2 + e2**2) 

    omega = phi_rdot + 2*e2*v_r*np.cos(e3_aux/2 - alpha) + (v_r**2*np.sin(e3) - e2*v_rdot)/fr + k2*np.sin(e3_aux/2)
    v = v_r*np.cos(e3) - omega*v_r*np.sin(e3_aux/2)/fr + k1*e1

    return np.array((v, omega))
    