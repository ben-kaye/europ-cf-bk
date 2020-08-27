import osqp
import numpy as np
import scipy as sp
from scipy import sparse

# x in Re(3x1) { p_x, p_y, phi }
# r in Re(6x1) { r_x, r_y, r_phi, r_phidot, r_v, r_vdot }
# u in Re(2x1) { v, omega }

sim_time = 10
step_size = 1e-4
Ts = 1e-3

k1 = 5
k2 = 5
max_turn = 1.5
v_min = 0.01
v_max = 12
delta = 0.4

x = np.array(0.3, 1.2, -np.radians(10)) # init state
r = np.array(2, 3, np.radians(60), 0.5, 1, 0) # init reference signal
p_o = np.array(-1, 6) # init object pos

# sim set up

N = sim_time/step_size
Ns = sim_time/Ts

xt = np.zeros(Ns)
rt = np.zeros(Ns)
ut = np.zeros(Ns)

# qp set up

ctrl_u = clf_controls(x, r, k1, k2)
ctrl_u = saturate_ctrls(u, max_turn, v_min, v_max)

Av = np.array((1, 0))
lv = v_min
uv = v_max

Aturn = np.array((0, 1))
lturn = -max_turn
uturn = max_turn

Abf = np.array((0, -1))
lbf = -np.inf
ubf = get_sig_bf_constr(x, ctrl_u, p_o, max_turn, delta)

A = sparse.vstack(Av, Aturn, Abf)
l = np.hstack(lv, lturn, lbf)
u = np.hstack(uv, uturn, ubf)

P = sparse.eye(2)
q = -2*P.dot(ctrl_u)

solver = osqp.OSQP()

solver.setup(P, q, A, l, u, warm_start = True, verbose = False)

for i in range(Ns):
    # solve
    res = solver.solve()

    if res.info.status != 'solved':
        raise ValueError('OSQP did not solve the problem!')

    ctrl_u = res.x

    # simulate
    for j in range(N/Ns):
        x = sim_step(x, ctrl_u, step_size)
        r = ref_step(r, step_size)

    # capture state
    xt[i] = x
    rt[i] = r
    ut[i] = ctrl_u

    # update
    ctrl_u = clf_controls(x, r, k1, k2)
    ctrl_u = saturate_ctrls(u, max_turn, v_min, v_max)

    u[2] = get_sig_bf_constr(x, ctrl_u, p_o, max_turn, delta)

    q = -2*P.dot(ctrl_u)

    solver.update(u=u, q=q)



def sim_step(x, u, step_sz):
    phi = x[2]

    Ax = np.array((np.cos(phi), 0),(np.sin(phi), 0),(0, 1))

    return x + step_sz*Ax.dot(u)

def ref_step(r, step_sz):
    phi_r = r[2]

    r_dot = np.array(np.cos(phi_r), 0),(np.sin(phi_r), 0),(0, 1))
    return r + step_sz*np.vstack(r_dot(r[:2]), 0, r[4], 0)  
    
def saturate_ctrls(u, max_turn, v_min, v_max):
    v = u[0]
    v = max(v, v_min)
    v = min(v, v_max)

    omega = u[1]
    omega = np.sign(omega)*max(abs(omega))

    return np.array(v, omega)

def get_sig_bf_constr(x, u, p_o, max_turn, delta):
    p = x[:1]
    phi = x[2]
    v = u[0]

    dir_v = np.array(cos(phi), sin(phi))

    p_xo = p_o - p

    cross_k = np.array((0, 1),(-1, 0)) 

    sign_r = np.sign(dir_v.dot(cross_k.dot(p_xo)))
    R_min = v/max_turn

    r = sign_r*R_min*cross_k.dot(dir_v)

    z = p_xo - r

    h = z.dot(z) - (delta + R_min)**2

    q = exp(1/h**2)

    sig = (q - 1)/(q + 1)

    # Abf = np.array((0, -1))
    ubf = sig*sign_r*max_turn

    return ubf

def clf_controls(x, r, k1, k2):
    p = x[:1]
    phi = x[2]
    p_r = r[:1]
    phi_r = r[2]
    phi_rdot = r[3]
    v_r = r[4]
    v_rdot = r[5]

    e3 = phi_r - phi
    c = np.cos(phi)
    s = np.sin(phi)
    e12 = np.array((c, s), (-s, c)).dot(p_r-p_x)

    alpha = np.atan(e12[1]/v_r)
    e3_aux = e3 + alpha

    fr = (v_r**2 + e12[1]**2) 

    omega = phi_rdot + 2*e12[1]*v_r*np.cos(e3_aux/2 - alpha) + (v_r**2*np.sin(e3) - e12[1]*v_rdot)/fr + k2*np.sin(e3_aux/2)
    v = v_r*np.cos(e3) - omega*v_r*np.sin(e3_aux/2)/fr + k1*e1

    return np.array(v, omega)
    