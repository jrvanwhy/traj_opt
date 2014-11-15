% Copyright (C) 2014 Johnathan Van Why
% See LICENSE.txt for details

addpath('..')

clear all
nodes = 50;
prob  = OptTool;

pos         = prob.newVar('pos', .5 * ones(nodes, 1));
T           = prob.newVar('T',   1);
[vel,accel] = diffTraj(pos, T);

prob.addCon(vel(1),   '==', 0);
prob.addCon(vel(end), '==', 0);
prob.addCon(pos(1),   '==', 0);
prob.addCon(pos(end), '==', 1);
prob.addCon(T,        '>=', 0);

% Do this variable last -- addCon's runtime scales with the number of variables.
u = prob.newVar('u', zeros(nodes, 1), -1, 1);
prob.addCon(accel, '==', u);

prob.addObj(T)

prob.solve

% Plot the optimal control input
plot(prob.getVar(u))
