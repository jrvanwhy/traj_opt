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

% This class represents a vectorized expression derived from symbolic
% inputs. It abstracts away differentiation and evaluation of
% expressions for the OptTool class.

classdef ExprEvaluator < handle
	methods
		% Constructor. This handles the necessary indexing
		% abstractions for OptTool.
		%
		% Parameters:
		%     expr   The Expr object that this will evaluate
		%     otool  A reference to the OptTool instance for symbolic
		%            variable access.
		function self = ExprEvaluator(expr, otool)
			% If the user passed in an anonymous function, switch out the variables
			% in the function for dummy variables.
			if ~isempty(expr.vars)
				fcn_vars = sym([otool.symvar_prefix 'var'], [numel(expr.vars) 1]);
				fcn_vars_cell = num2cell(fcn_vars);
				expr.expr = expr.expr(fcn_vars_cell{:});
			else
				fcn_vars = otool.vars;
			end

			% Generate the anonymous function for evaluating the objective itself
			self.fcn = matlabFunction(sym(expr.expr), 'vars', {fcn_vars});

			% Create a sparse representation of the jacobian of each term
			% in this vectorized expression. i is in terms of the term's subexpressions,
			% j is in terms of fcn_vars, and s is per unique derivative value
			[jac_i, jac_j, jac_s] = find(jacobian(expr.expr, fcn_vars));

			% Indices used in this vectorized expression
			if isempty(expr.idxs)
				expr.idxs = 1:max(otool.var_sizes(jac_j));
			end

			% Generate the variable map which maps optimization variables into inputs for
			% self.fcn and self.jacfcn
			if ~isempty(expr.vars)
				self.var_map = [];

				% Handle the variables one at a time. We'll expand var_map as necessary
				% while running this loop.
				for var_idxs = expr.idxs
					% Handle a for loop peculiarity with cell arrays
					var_idxs = var_idxs{1};

					% Expand self.var_map if necessary
					if numel(var_idxs) > size(self.var_map, 2)
						self.var_map = repmat(self.var_map, 1, numel(var_idxs));
					end

					% Append to self.var_map
					self.var_map = [self.var_map; var_idxs];
				end
			else
				self.var_map = bsxfun(@plus, otool.var_start_idxs, (otool.var_sizes > 1).' * (expr.idxs - 1));
			end

			% Compute and store the number of outputs for the user
			self.num_outs = size(self.var_map, 2);

			% Generate the jacobian function.
			self.jacfcn = matlabFunction(jac_s(:), 'vars', {fcn_vars});

			% Generate the indexing lists for the sparse representation of the jacobian
			full_jac_i = bsxfun(@plus, jac_i(:) - 1, numel(expr.expr) * (1:size(self.var_map, 2)));
			self.jac_i = full_jac_i(:);
			full_jac_j = self.var_map(jac_j, :);
			self.jac_j = full_jac_j(:);
		end

		% Evaluator; evaluates the function at the given point
		function val = eval(self, x)
			val = self.fcn(x(self.var_map));
			val = val(:);
		end

		% Jacobian evaluator; evaluates the jacobian of the function
		% at the given point.
		function jac = eval_jac(self, x)
			jac_s = self.jacfcn(x(self.var_map));
			jac = sparse(self.jac_i, self.jac_j, jac_s(:), max([0; self.jac_i]), numel(x));
		end
	end

	properties
		% The number of outputs for this vectorized expression
		num_outs@double

		% Map between the optimization problem's variables and the opttool variables
		var_map@double

		% Anonymous function for evaluating this vectorized expression
		fcn

		% Anonymous function for evaluating the jacobian entries
		% as well as the coordinates of the entries
		jacfcn
		jac_i@double
		jac_j@double
	end
end
