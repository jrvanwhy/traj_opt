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

	% Generate the fields for the function values themselves
	vecfcn.params = param_nums;
	vecfcn.fcn    = matlabFunction(expr, 'vars', {sym_params});

	% Exit now if we don't care about derivatives
	if ~derivs
		return
	end

	% Generate the fields for the jacobian
	% If sym_params is empty, work around a crash in jacobian()
	if isempty(sym_params)
		jac_expr = sym(zeros(1, 0));
	else
		jac_expr = jacobian(expr, sym_params);
	end
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
