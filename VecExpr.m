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

classdef VecExpr < handle
	methods
		% Constructor. This handles the necessary indexing
		% abstractions for OptTool.
		%
		% Parameters:
		%     otool  A reference to the OptTool instance for symbolic
		%            variable access.
		%     expr   Symbolic expression, as in OptTool's functions
		%     idxs   Indices list, as in OptTool's interface
		%     vars   Variables list, also as in OptTool's API
		function self = VecExpr(otool, expr, idxs, vars)
			% TODO: Handle idxs and vars-based inputs...

			% Basic case -- just a symbolic expression with
			% no indexing
			self.fcn = matlabFunction(expr, 'vars', {otool.vars});
		end
	end

	properties
		% Anonymous function for evaluating this vectorized expression
		fcn
	end
end
