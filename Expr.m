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

% This class represents a vectorized expression in an optimization.
% The expression may be a simple symbolic expression (operating on symbolic variable
% groups), a symbolic expression with indexing into its output, or an expression
% defined by an anonymous function with different variable group indexing for its
% different inputs.

classdef Expr < handle
	methods
		% Constructor. Expr encapsulates everything required for to represent
		% a symbolic vectorized function.
		%
		% Parameters:
		%     expr   Symbolic expression operating on symbolic variable groups or
		%            an anonymous function
		%     idxs   Indices list (optional) -- if provided, gives indices
		%            into the expression group's output or (if vars is provided)
		%            indices for the inputs into the anonymous function expr
		%     vars   Variables list (optional) -- inputs to expr if expr
		%            is a function handle
		function self = Expr(expr, idxs, vars)
			if nargin < 3
				vars = [];

				if nargin < 2
					idxs = [];
				end
			end

			self.expr = expr;
			self.idxs = idxs;
			self.vars = vars;
		end
	end

	properties
		% The values passed into this structure at construction time.
		% These are empty if the relevant variable was not passed to the constructor
		expr
		idxs
		vars
	end
end
