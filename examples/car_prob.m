% Copyright (C) 2016 Johnathan Van Why
% See LICENSE.txt for details

addpath('..')

clear all

nlp = OptTool;

t = nlp.newVar('t', 1);

nodes = 10;
dt = t / (nodes - 1);

pos = nlp.newVar('pos', [0 linspace(0, 1, nodes) 1]);
acc = nlp.newVar('acc', zeros(nodes, 1));

nlp.addCon(pos, '==', 0, 2);
nlp.addCon(pos, '==', 1, nodes+1);
nlp.addCon(acc, '>=', -1);
nlp.addCon(acc, '<=', 1);

nlp.addCon(Expr(@(xi, xc, xf, u, t) xf + xi - 2*xc - (t/(nodes-1))^2 * u, {1:nodes, 2:nodes+1, 3:nodes+2, 1:nodes 1}, [pos pos pos acc t]), '==', 0);

nlp.addCon(Expr(@(xi, xf) xf - xi, {1 3}, [pos pos]), '==', 0);
nlp.addCon(Expr(@(xi, xf) xf - xi, {nodes nodes+2}, [pos pos]), '==', 0);

nlp.addObj(t);

nlp.solve;
