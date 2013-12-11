% This function evaluates the given vecfcn with the given parameters

function vals = opt_eval_vecfcn(vecfcn, params)
	% Create the parameters matrix with the correct dimensions
	vecparams    = zeros(size(vecfcn.params));

	% Make sure we don't crash if params is symbolic
	if ~isnumeric(params)
		vecparams = sym(vecparams);
	end

	vecparams(:) = params(vecfcn.params);

	% Call the function (it's that easy!)
	vals = vecfcn.fcn(vecparams);
end
