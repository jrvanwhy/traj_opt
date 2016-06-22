% Copyright (C) 2016 Johnathan Van Why
% See LICENSE.txt for details

addpath('..')

clear all

nlp = OptTool;

x = nlp.newVar('x', exp(1));
y = nlp.newVar('y', pi);

nlp.addObj([(1 - x)^2, 100 * (y - x^2)^2]);

nlp.solve;
