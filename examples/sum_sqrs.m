% Copyright (C) 2016 Johnathan Van Why
% See LICENSE.txt for details

addpath('..')

clear all

nlp = OptTool;

x = nlp.newVar('x', 1:10);
y = nlp.newVar('y', 1);

nlp.addObj((x - y)^2, 1:5);
nlp.addObj(x^2, 6:10);
nlp.addObj((y - pi)^2);

%nlp.addObj((x - y)^2);

nlp.solve;
