% Copyright (C) 2016 Johnathan Van Why
% See LICENSE.txt for details

addpath('..')

clear all

nlp = OptTool;

x_size = 10;
x = nlp.newVar('x', zeros(x_size, 1));

nlp.addObj(@(x1, x2) 100 * (x2 - x1^2)^2 + (x1 - 1)^2, {1:x_size-1, 2:x_size}, [x x]);

nlp.addCon(x.^2, '<=', 1, 1:5);

nlp.solve;
