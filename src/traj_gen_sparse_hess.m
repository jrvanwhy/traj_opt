% This function creates a sparse hessian. You need to give it the parameters
% the hessian is in terms of, the lambdas (may be symbolic or numeric), and
% the jacobian of the value

function hess = traj_gen_sparse_hess(params, lambdas, jac)
	% Multiply the lambdas and the jacobian to get the jacobian of the
	% sum of the lambdas times the functions
	sum_jac = lambdas.' * jac

	% Extract the nonzero column numbers and values
	[~,nonzero_cols,nonzero_elems] = find(sum_jac)

	% Compute the hessian (only nonzero rows)
	red_hess = jacobian(nonzero_elems, params)

	% Convert reduced hessian to a sparse hessian
	hess.mn = numel(params)
	[hess.i,hess.j,hess.s] = find(red_hess)

	% Create non-reduced hessian by filling in the zero rows
	hess.i = nonzero_cols(hess.i)
end
