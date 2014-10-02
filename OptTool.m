% This is a simple tool for formulating and solving nonlinear optimization problems.
% It does a few things to try to increase performance, but won't be terribly fast.

classdef OptTool < handle
	methods
		% Adds a variable.
		%
		% params
		%     name    The variable's name. Must be a valid symbolic variable.
		%     initVal The initial value for the variable; the variable will have the same dimensions as this vector.
		%     lb      The variable's lower bound (optional; use -inf to specify no lower bound for a particular element).
		%     ub      Upper bound. Optional; use inf to specify no upper bound for a particular element of this variable.
		%
		% Returns a symbolic expression for the variable
		function expr = newVar(this, name, initVal, lb, ub)
			% There's not much error checking here, but this one is particularly easy
			if nargin < 3
				error('Not enough arguments')
			end

			% Basic diagnostics for the user
			disp(['Creating variable ' name ' of size ' num2str(numel(initVal))])

			% Create the variable. If it's scalar, create it as a scalar (this leads to nicer expressions)
			if numel(initVal) == 1
				expr = sym(name, 'real');
			else
				expr = sym(name, [numel(initVal) 1]);
			end

			% Convenience variable -- all ones column vector of the same size as this variable
			szones = ones(numel(initVal), 1);

			% Copy over the lower and upper bounds, if supplied -- if not, generate them.
			% Watch out for scalars!
			if nargin < 4 || isempty(lb)
				lb = -inf * szones;
			elseif numel(lb) == 1
				lb = lb * szones;
			end

			if nargin < 5 || isempty(ub)
				ub = inf * szones;
			elseif numel(ub) == 1
				ub = ub * szones;
			end

			this.lb = [this.lb; lb];
			this.ub = [this.ub; ub];

			% Append to variable-related members
			this.varmap(1:2, end+1) = [expr(1) 1+numel(this.vars)];
			this.vars               = [this.vars; expr   ];
			this.x0                 = [this.x0;   initVal];
		end

		% Adds an expression to the problem objective.
		%
		% params:
		%     expr The expression to add to the objective
		function addObj(this, expr)
			% User diagnostics
			disp('Adding objective')

			% Just add it to the class's member. Unlike the constraints, we don't
			% calculate the gradient yet. With the objective, adding new terms does not
			% affect gradient size, so it is most efficiently done in one big sweep.
			this.obj = this.obj + expr;
		end

		% Add is a constraint. Type should be one of '<=', '==', or '>='
		%
		% params:
		%     lexpr The expression on the left side of the constraint
		%     type  The type of constraint
		%     rexpr The expression on the right side of the constraint
		function addCon(this, lexpr, type, rexpr)
			disp(['Adding constraint of type ' type])

			% Convert the constraint to the form g(x) <= 0 or h(x) == 0
			switch type
				case '<='
					expr = lexpr - rexpr;

				case '=='
					expr = lexpr - rexpr;

				case '>='
					expr = rexpr - lexpr;

				otherwise
					error(['Constraint type ''' type ''' invalid'])
			end

			% Append this constraint to the relevant expression.
			% Also grab the "current" jacobian based on the constraint type
			if type == '=='
				this.ceq = [this.ceq; expr];
				curjac   = this.ceqjac;
			else
				this.c   = [this.c;   expr];
				curjac   = this.cjac;
			end

			% Begin Jacobian calculation -- evaluate the sparse jacobian for this constraint only first
			spjac = this.sparse_jacobian(expr, this.vars);

			% Shift over the new jacobian elements to lie "past" the current elements (shifting to the right horizontally).
			spjac.j = spjac.j + curjac.n;

			% Append/update jacobian expressions
			curjac.i = [curjac.i; spjac.i(:)];
			curjac.j = [curjac.j; spjac.j(:)];
			curjac.s = [curjac.s; spjac.s(:)];
			curjac.n = curjac.n + spjac.n;

			% Apply the updated jacobian to "this"
			if type == '=='
				this.ceqjac = curjac;
			else
				this.cjac = curjac;
			end
		end

		% Grabs the solution for the given variable. Must be an original variable returned by newVar()
		%
		% params:
		%     var The symbolic variable to look for
		%
		% returns:
		%     soln: The solution for this variable
		function soln = getVar(this, var)
			soln = this.soln(this.varmap(2,find(logical(this.varmap(1,:) == var(1)))) + (1:numel(var)) - 1);
		end

		% Evaluates a sparse function
		function val = eval_sparse(this, spfcn, vars)
			% Evaluate the elements
			spfcn.s = spfcn.s(vars);

			% Create a sparse matrix representing the result
			val = sparse(spfcn.i, spfcn.j, spfcn.s, spfcn.m, spfcn.n);
		end

		% Computes a sparse symbolic jacobian (transposed, because fmincon wants gradients transposed)
		function spjac = sparse_jacobian(this, expr, vars)
			disp('Evaluating full jacobian')
			jac = jacobian(expr, vars);

			disp('Locating nonzero jacobian entries')
			[spjac.j,spjac.i,spjac.s] = find(jac);
			disp([num2str(numel(spjac.i)) ' entries found'])

			disp('Storing jacobian dimensions')
			spjac.n = size(jac, 1);
			spjac.m = size(jac, 2);
			disp(['Dimensions: ' num2str(spjac.n) ' by ' num2str(spjac.m)])
		end

		% Lets the user alter the Fmincon options
		function setOptions(this, varargin)
			this.options = optimset(this.options, varargin{:});
		end

		% Objective function called by fmincon
		function [obj,gobj] = fmincon_obj(this, x)
			% Simply forward the obj calculation to the objfcn handle
			obj = this.objfcn(x);

			% Do a sparse gradient evaluation, if necessary
			if nargout > 1
				gobj = this.eval_sparse(this.objjac, x);
			end
		end

		% Nonlinear constraints function called by fmincon
		function [c,ceq,gc,gceq] = fmincon_cons(this, x)
			% Simply forward constraint evaluations to the handles (they're dense)
			c = this.cfcn(x);

			% Quit early if possible
			if nargout <= 1
				return
			end

			% Same for ceq
			ceq = this.ceqfcn(x);
			if nargout <= 2
				return
			end

			% Evaluate sparse gradients similar to the objective
			gc = this.eval_sparse(this.cjac, x);
			if nargout <= 3
				return
			end

			gceq = this.eval_sparse(this.ceqjac, x);
		end

		% Solves the optimization problem.
		% The result should be retrieved using the getSoln() function
		function solve(this)
			disp('Beginning opttool.solve()')

			% Generate the objective function
			disp('Generating objective function')
			this.objfcn = matlabFunction(this.obj, 'vars', {this.vars});

			% Generate the jacobian of the objective function
			disp('Generating objective jacobian')
			this.objjac = this.sparse_jacobian(this.obj, this.vars);

			disp('Generating anonymous jacobian evaluation function')
			this.objjac.s = matlabFunction(this.objjac.s, 'vars', {this.vars});

			% Generate similar functions for the constraint gradients
			disp('Generating inequality constraint function')
			this.cfcn   = matlabFunction(this.c,   'vars', {this.vars});

			disp('Generating equality constraint function')
			this.ceqfcn = matlabFunction(this.ceq, 'vars', {this.vars});

			disp('Finalizing inequality constraint jacobian')
			this.cjac.s = matlabFunction(this.cjac.s, 'vars', {this.vars});
			this.cjac.m = numel(this.vars);

			disp('Finalizing equality constraint jacobian')
			this.ceqjac.s = matlabFunction(this.ceqjac.s, 'vars', {this.vars});
			this.ceqjac.m = numel(this.vars);

			% Finally, we call fmincon
			disp('Starting Fmincon')
			this.soln = fmincon(@this.fmincon_obj, this.x0, [], [], [], [], this.lb, this.ub, @this.fmincon_cons, this.options);
		end
	end % methods

	properties
		% Lower and upper bounds for the decision variables
		lb@double
		ub@double

		% The objective expression (symbolic) and function (an anonymous function) as well as other related values
		objfcn
		objjac
		obj@sym = sym(0)

		% Constraint and constraint jacobian expressions (like the objective ones)
		c@sym   = sym([])
		ceq@sym = sym([])
		cfcn
		ceqfcn
		cjac    = struct('i', [], 'j', [], 's', sym([]), 'm', 0, 'n', 0);
		ceqjac  = struct('i', [], 'j', [], 's', sym([]), 'm', 0, 'n', 0);

		% Fmincon options
		options = optimset('Algorithm',  'interior-point', ...
		                   'Display',    'iter',           ...
		                   'GradObj',    'on',             ...
		                   'GradConstr', 'on');

		% The solution fmincon returns
		soln

		% The objective decision variables, as a symbolic column vector
		vars@sym = sym([])

		% This essentially forms a dictionary, matching the first element of variables with their starting index.
		varmap@sym = sym([])

		% Initial guess for the problem
		x0@double
	end % properties
end % OptTool
