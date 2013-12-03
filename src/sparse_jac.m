% This function computes the sparse jacobian of the given expression
% In other words, it is a sparse version of jacobian() (part of the symbolic toolkit)
% expr must be a column vector (or scalar) symbolic expression, and params must be a symbolic column vector of variables.
% It returns a structure containing all the arguments to pass to sparse() to generate
% the numerical symbolic matrix (the "s" function being symbolic, though).

% Full is used internally and should not be user-specified

function jac = sparse_jac(expr, params, non_params, full)
	% exprvars is the list of optimization parameters in this 
	if nargin < 4
		% Optimization for the first call
		exprvars = params;
	else
		% The fastest way I've found to do this
		exprvars = symvar(subs(symvar(expr), non_params, zeros(size(non_params)))).';
	end

	% Number of values in expr
	n_outs = numel(expr);

	% Check if this is too big to compute directly; if so, split it in two
	% This value should be tuned better; currently, it's only a rough estimate
	if n_outs > 1 && n_outs * numel(exprvars) > 10000
		% The middle output
		mid_idx = floor(n_outs/2);
		jac = sparse_vcat(sparse_jac(expr(    1:mid_idx), params, non_params, 1), ...
		                  sparse_jac(expr(mid_idx+1:end), params, non_params, 1));

		% Quit early -- we have our answer already
		return
	end

	% Okay, it's small enough to compute in full

	% Set the size of jac
	jac.m = n_outs;
	jac.n = numel(params);

	% Compute the location of each parameter in params
	varlocs = locate_sym(exprvars, params);

	% Do the jacobian
	jac_full = jacobian(expr, exprvars);

	% Find the nonzero values
	[jac.i,cols,jac.s] = find(jac_full);

	% Correct the columns to account for the fact that we did a smaller than necessary jacobian
	jac.j = varlocs(cols);
end

% This function vertically concatenates two sparse matrices
function matout = sparse_vcat(mat1, mat2)
	% Set matout size
	matout.m = mat1.m + mat2.m;
	matout.n = mat1.n;

	% Append the values, adjusting row numbers
	matout.i = [mat1.i; mat2.i + mat1.m];
	matout.j = [mat1.j; mat2.j];
	matout.s = [mat1.s; mat2.s];
end

% This function locates the variables in A in B, and returns their locations.
% It is equivalent to [~,varlocs] = ismember(A, B)
% but also works for symbolic expressions (ismember doesn't work in 2012b, this is probably a bug)
function varlocs = locate_sym(A, B)
	% Preallocation
	varlocs = zeros(size(A));

	% Iterate through one-by-one
	parfor iterA = 1:numel(A)
		varlocs(iterA) = find(B == A(iterA));
	end
end
