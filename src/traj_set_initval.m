% This function gives a parameter or parameters an initial
% value or initial values
% The number of initial values must be 1 or equal to the number of parameters

function strout = traj_set_initval(params, values)
	% Error checking
	if numel(values) ~= 1 && numel(values) ~= numel(params)
		error('Number of values must be 1 or equal to the number of parameters')
	end

	% Initialize output structure
	strout.version     = traj_version();
	strout.struct_type = 'initval';

	% Copy over parameters and values (replicating values if necessary)
	strout.params                  = params(:);
	strout.values(1:numel(params)) = values(:);

	% Reorient strout.values to be a column vector
	strout.values = strout.values(:);
end
