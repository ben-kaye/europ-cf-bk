import osqp
import numpy as np
import scipy as sp
from scipy import sparse
import qpfuncs as qp

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

x = np.array((0.3, 1.2, -np.radians(10))) # init state
r = np.array((2, 3, np.radians(60), 0.5, 1, 0)) # init reference signal
p_o = np.array((-1, 6)) # init object pos

# sim set up

N = int(sim_time/step_size)
Ns = int(sim_time/Ts)

xt = np.zeros(Ns)
rt = np.zeros(Ns)
ut = np.zeros(Ns)

# qp set up

ctrl_u = qp.clf_controls(x, r, k1, k2)
ctrl_u = qp.saturate_ctrls(ctrl_u, max_turn, v_min, v_max)

Av = np.array((1, 0))
lv = v_min
uv = v_max

Aturn = np.array((0, 1))
lturn = -max_turn
uturn = max_turn

Abf = sparse.csr_matrix(np.array((0, -1)))
lbf = -np.inf
ubf = qp.get_sig_bf_constr(x, ctrl_u, p_o, max_turn, delta)

A = sparse.vstack([Av, Aturn, Abf], format='csc')
l = np.hstack([lv, lturn, lbf])
u = np.hstack([uv, uturn, ubf])

P = sparse.eye(2, format="csc")
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
    for j in range(int(N/Ns)):
        x = qp.sim_step(x, ctrl_u, step_size)
        r = qp.ref_step(r, step_size)

    # capture state
    xt[i] = x
    rt[i] = r
    ut[i] = ctrl_u

    # update
    ctrl_u = qp.clf_controls(x, r, k1, k2)
    ctrl_u = qp.saturate_ctrls(u, max_turn, v_min, v_max)

    u[2] = qp.get_sig_bf_constr(x, ctrl_u, p_o, max_turn, delta)

    q = -2*P.dot(ctrl_u)

    solver.update(u=u, q=q)