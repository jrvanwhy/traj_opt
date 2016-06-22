% Copyright (C) 2016 Johnathan Van Why
% See LICENSE.txt for details

addpath('..')

clear all

nlp = OptTool;

x = nlp.newVar('x', 0);
y = nlp.newVar('y', 0);

nlp.addObj((1 - x)^2 + 10 * (y - x^2)^2);

nlp.solve;
