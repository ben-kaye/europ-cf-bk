# EUROP Project: Investigating optimal control of quad-copters

Notes: constr_ca_dl.m broken and not sure why?

## Drone Collision Avoidance
unconstrained_mpc.m &nbsp; &nbsp; &nbsp; MPC by solving ARE<br/>
constrained_mpc.m &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; MPC with speed constraint<br/>
constrained_mpc_ca.m &nbsp; &nbsp; &nbsp; &nbsp; MPC with linearised CA constraint<br/>
robust_mpc_ca.m &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp;Improved MPC with linearised CA constraint<br/>
ndrone_mpcca.m &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; N-Drones MPC CA<br/>

## Simulink
quadcop_simulator.slx <br/>
crazyflie_controller_ak.mdl
