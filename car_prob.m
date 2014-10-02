clear all
nodes = 100;
prob  = OptTool;

pos_var = prob.newVar('pos', .5 * ones(nodes-4, 1));
pos = [0; 0; pos_var; 1; 1];
t   = prob.newVar('t', 1);

dt = t/(nodes-1);

vel = (pos(2:end) - pos(1:end-1))/dt;

accel = (vel(2:end) - vel(1:end-1))/dt;

prob.addObj(t)

prob.addCon(accel, '<=',  1)
prob.addCon(accel, '>=', -1)

prob.solve

pos_soln = prob.getVar(pos_var);
