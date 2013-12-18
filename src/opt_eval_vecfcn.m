% This function evaluates the given vecfcn with the given parameters

function vals = opt_eval_vecfcn(vecfcn, params)
	% Create the parameters matrix with the correct dimensions
	vecparams    = zeros(size(vecfcn.param_nums));

	% Initialize vals to the constant values
	vals = vecfcn.fcn.const_outs;

	% Make sure we don't crash if params is symbolic
	if ~isnumeric(params)
		vecparams = sym(vecparams);
		vals      = sym(vals);
	end

	vecparams(:) = params(vecfcn.param_nums);

	% Call the function (it's that easy!) and write its outputs
	% into the correct rows in vals
	if numel(vecparams) > 0
		vals(vecfcn.fcn.nc_outs,:) = vecfcn.fcn.nc_fcn(vecparams);
	end
end
