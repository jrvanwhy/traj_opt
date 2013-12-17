% This function creates a new "vectorized function"
% A vectorized function is just a symbolic function that gets applied
% to a large number of parameters at a time. This library deals
% with calling the function appropriately and determining
% its derivatives
% The derivs argument controls whether derivatives should be computed
% It defaults to true if not specified

function vecfcn = opt_create_vecfcn(expr, param_nums, sym_params, derivs)
	% Default derivs to true if not given
	if nargin < 4
		derivs = true;
	end

	% Initialize the structure
	vecfcn.n_vals = numel(expr);

	% Handle the jacobian first, as it is used for simplification
	if derivs
		% Take the jacobian, working around a crash in jacobian()
		% that occurs when the parameters list is empty
		if isempty(sym_params)
			jac_expr = sym(zeros(vecfcn.n_vals, 0));
		else
			jac_expr = jacobian(expr, sym_params);
		end

		% Identify parameters used in expr
		% To do this, we find nonzero columns in the jacobian
		% First, we identify nonzero cells in the jacobian
		nz_cells = logical(jac_expr ~= 0);

		% Then sum each column and compare with 0 to find nonzero columns
		nz_cols = (sum(nz_cells, 1)) ~= 0;

		% Reduce the parameters list to only include those parameters actually in
		% the function
		param_nums = param_nums(nz_cols,:);
		sym_params = sym_params(nz_cols);
	end

	% Generate the fields for the function values themselves
	vecfcn.params = param_nums;
	vecfcn.fcn    = matlabFunction(expr, 'vars', {sym_params});

	% Exit now if we don't care about derivatives
	if ~derivs
		return
	end

	% The Jacobian was generated earlier, but for the full list of
	% symbolic parameters. To get the jacobian function for the
	% vecfcn structure, we need to remove the columns that don't correspond
	% with any variables
	jac_expr = jac_expr(:, nz_cols);

	% Create the jacobian anonymous function
	vecfcn.jac_fcn = matlabFunction(jac_expr, 'vars', {sym_params});

	% Generate the hessian. This requires generating the relevant symbolic lambda values
	lambdas = sym('h', [numel(expr) 1]);

	% Correctly multiply the lambdas with the jacobian to get the hessian of the lagrangian
	% Again, we need to work around symbolic engine issues if lambdas is empty
	if isempty(lambdas)
		hess_expr = sym(zeros(0, 0));
	else
		hess_expr = jacobian(lambdas.' * jac_expr, sym_params);
	end
	vecfcn.hess_fcn = matlabFunction(hess_expr, 'vars', {sym_params, lambdas});

	% Clean up after ourself
	syms lambdas clear
end
