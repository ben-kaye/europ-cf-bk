# EUROP Project: Investigating optimal control of quad-copters

Notes: constr_ca_dl.m broken and not sure why?

## Drone Collision Avoidance
unconstrained_mpc.m &emsp; &emsp; MPC by solving ARE<br/>
constrained_mpc.m &emsp; &emsp; MPC with speed constraint<br/>
constrained_mpc_ca.m &emsp; &ensp; MPC with linearised CA constraint<br/>
robust_mpc_ca.m &emsp; &emsp; &ensp; Improved MPC with linearised CA constraint<br/>
ndrone_mpcca.m &emsp; &emsp; N-Drones MPC CA<br/>

## Simulink
quadcop_simulator.slx <br/>
crazyflie_controller_ak.mdl
