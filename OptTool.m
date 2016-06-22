% Copyright (C) 2016 Johnathan Van Why
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2.1
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

% The license may be found in LICENSE.txt

% This is a simple tool for formulating and solving nonlinear optimization problems.
% It attempts to achieve reasonably high performance without an excessive amount of
% complexity.

classdef OptTool < handle
	methods
		% Constructor; initializes this OptTool instance with the given
		% settings (if any are provided)
		%
		% Parameters:
		%     symvar_prefix  Prefix for the symbolic variables.
		%                    Defaults to 'opt_'. The prefix is used
		%                    to prevent name conflicts (especially
		%                    with MuPAD reserved names).
		function self = OptTool(symvar_prefix)
			% Default argument handling
			if nargin < 1
				symvar_prefix = 'opt_';
			end

			% Set the provided properties
			self.symvar_prefix = symvar_prefix;
		end

		% Adds a variable.
		%
		% Parameters:
		%     name     The variable's name (affects the generated symbolic expression)
		%     initVal  Initial value fonr the variable. The variable will have the same dimensions as this vector.
		%     lb       The variable's lower bound (optional -- use -inf to specify no lower bound for a particular element)
		%     ub       The variable's upper bound (optional -- use inf to specify no upper bound for a particular element)
		%
		% Returns a scalar symbolic expression representing this variable.
		function expr = newVar(self, name, initVal, lb, ub)
			% Basic diagnostics for the user
			disp(['Creating variable ' name ' of size ' num2str(numel(initVal))]);

			% Create the symbolic variable
			expr = sym([self.symvar_prefix name], 'real');

			% Append this variable to the variable-related members
			if ~isempty(self.vars)
				self.var_start_idxs(end+1, 1) = self.var_start_idxs(end) + self.var_sizes(end);
			else
				self.var_start_idxs = 1;
			end
			self.var_sizes(end+1) = numel(initVal);
			self.vars(end+1, 1) = expr;
			self.x0 = [self.x0; initVal(:)];
		end

		% Adds an expression to the problem objective.
		% If this expression is non-scalar, then its sum will be added
		% to the problem objective.
		%
		% Parameters:
		%     expr  Symbolic expression for this term or an anonymous function/
		%           function handle representing this expression.
		%     idxs  Indices into the expression to include in the sum, or
		%           (if vars is also provided) index lists for each variable
		%     vars  A cell array of variables corresponding with idxs or with
		%           the parameters of expr
		function addObj(self, varargin)
			% Diagnostic for the user
			disp('Adding objective');

			% Behind the scenes, make the VecExpr class do all the hard work
			self.objs(end+1) = VecExpr(self, varargin{:});
		end

		% General VecExpr array evaluator; for internal use.
		%
		% Parameters:
		%     fcns      A cell array of VecExpr to evaluate
		%     x         The point to evaluate at
		%     calc_jac  Whether or not to compute the jacobian. If false,
		%               will return an empty jacobian. This parameter is a
		%               convenience for the call sites
		%
		% Returns:
		%     vals  All of the function outputs, vertically concatenated
		%     jac   The jacobian of the functions in terms of x
		function [vals,jac] = evalFcns(self, fcns, x, calc_jac)
			vals = [];
			jac = sparse(0, numel(x));

			for fcn = fcns
				vals = [vals; fcn.eval(x)];

				if calc_jac
					jac = [jac; fcn.eval_jac(x)];
				end
			end
		end

		% Objective function for the fmincon solver. This function satisfies
		% the interface required by fmincon
		function [obj,gobj] = fminconObj(self, x)
			% Compute all of the objective values
			[obj,jobj] = self.evalFcns(self.objs, x, nargout >= 2);

			% Sum up the objective values
			obj = sum(obj);
			if nargout >= 2
				gobj = sum(jobj, 1).';
			end
		end

		% Set Fmincon options. This accepts the same arguments as optimset
		function setOptions(self, varargin)
			self.options = optimset(self.options, varargin{:});
		end

		% Solves the nonlinear optimization problem. This needs to be called
		% after problem setup (adding variables, objectives, constraints, etc...)
		% and before querying solution values.
		function solve(self)
			% User diagnostic
			disp('Beginning OptTool.solve()');

			disp('Starting Fmincon');
			self.soln = fmincon(@self.fminconObj, self.x0, [], [], [], [], [], [], [], self.options);
		end
	end

	properties
		% Prefix for the symbolic expressions.
		symvar_prefix@char

		% All of the symbolic variables associated with this problem instance.
		vars@sym = sym([])

		% Starting indices and sizes for the variables
		var_start_idxs = []
		var_sizes = []

		% Initial guess for the solver for this problem
		x0@double = []

		% List of objectives (which will be summed to evaluate the final
		% objective)
		objs@VecExpr

		% Fmincon options structure
		options = optimset('Algorithm',           'interior-point', ...
		                   'DerivativeCheck',     'off',            ...
		                   'Display',             'iter',           ...
		                   'GradObj',             'on',             ...
		                   'SubproblemAlgorithm', 'cg');

		% The solution value (as returned by the nonlinear programming solver)
		soln@double
	end
end
