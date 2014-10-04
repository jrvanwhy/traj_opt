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

% Setting the Hessian option to user-supplied uses a Hessian*vector
% multiply function inside OptTool. For some reason, this will cause
% fmincon to fail on some problems (compared to the default setting,
% 'fin-diff-grads'), but usually improves fmincon's performance.
% Therefore, it defaults to off, but is a setting worth considering.
% For this problem, however, this has a noticeable benefit
prob.setOptions('Hessian', 'user-supplied')

prob.solve
