% Copyright (C) 2014 Johnathan Van Why
% See LICENSE.txt for details

addpath('..')

clear all

prob = OptTool;

x = prob.newVar('x', 0);
y = prob.newVar('y', 0);

prob.addObj((1 - x)^2 + 100 * (y - x^2)^2)

prob.solve
