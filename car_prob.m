clear all
nodes = 1000;
prob  = OptTool;

pos         = prob.newVar('pos', .5 * ones(nodes, 1));
u           = prob.newVar('u',   zeros(nodes, 1), -1, 1);
T           = prob.newVar('T',   1                  );
[vel,accel] = diffTraj(pos, T);

prob.addCon(accel,    '==', u);
prob.addCon(vel(1),   '==', 0);
prob.addCon(vel(end), '==', 0);
prob.addCon(pos(1),   '==', 0);
prob.addCon(pos(end), '==', 1);
prob.addCon(T,        '>=', 0);

prob.addObj(T)

prob.solve

pos_soln = prob.getVar(pos);
