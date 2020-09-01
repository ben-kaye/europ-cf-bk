import osqp
import numpy as np
import scipy as sp
from scipy import sparse
from scipy.io import savemat
import qpfuncs as qp

# x in Re(3x1) { p_x, p_y, phi }
# r in Re(6x1) { r_x, r_y, r_phi, r_phidot, r_v, r_vdot } {}
# u in Re(2x1) { v, omega }

# Implementing min-norm QP with sigmoid constraint based BF

sim_time = 10
step_size = 1e-4
Ts = 1e-3

k1 = 3
k2 = 5
max_turn = 1.5
v_min = 0.01
v_max = 12
delta = 0.4
a = 1


x = np.array((3, 2, -np.radians(30))) # init state
r = np.array((2, 3, np.radians(60), 0.5, 1, 0)) # init reference signal
p_o = np.array((-1, 6)) # init object pos

# sim set up

N = int(sim_time/step_size)
Ns = int(sim_time/Ts)

xt = np.zeros([Ns,3])
rt = np.zeros([Ns,3])
ut = np.zeros([Ns,2])

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
ubf = qp.get_sig_bf_constr(x, ctrl_u, p_o, max_turn, delta,a)

A = sparse.vstack([Av, Aturn, Abf], format='csc')
l = np.hstack([lv, lturn, lbf])
u = np.hstack([uv, uturn, ubf])

P = sparse.eye(2, format="csc")
q = -2*P.dot(ctrl_u)

solver = osqp.OSQP()

solver.setup(P, q, A, l, u, warm_start = False, verbose = False)

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
    rt[i] = r[:3]
    ut[i] = ctrl_u

    ctrl_temp = ctrl_u

    # update
    ctrl_u = qp.clf_controls(x, r, k1, k2)
    ctrl_u = qp.saturate_ctrls(ctrl_u, max_turn, v_min, v_max)

    u[2] = qp.get_sig_bf_constr(x, ctrl_temp, p_o, max_turn, delta, a) # note using v_temp instead of ctrl_u

    q = -2*P.dot(ctrl_u)

    solver.update(u=u, q=q)

# write to .csv for matlab plotting
sp.io.savemat('H:\\Files\\EUROP-MATLAB\\Python-CBF\\result.mat', {'x_t': xt, 'r_t': rt, 'u_t': ut, 'p_o': p_o, 'delta': delta})

print('Solved')