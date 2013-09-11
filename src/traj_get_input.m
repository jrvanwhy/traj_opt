% This function returns the input at a given time. Subinputs specifies which part of the input to return (if given).
% Scenario must have been run through traj_gen_numerical prior to calling this function

function input = traj_get_input(scenario, t, subinputs)
	% If subinputs was not specified, return the whole input
	if nargin < 3
		subinputs = 1:size(scenario.inputs.inputs, 1);
	end

	% Time values for all input points
	x = linspace(0, scenario.num_duration, size(scenario.inputs.inputs, 2));

	% All the inputs (just the subinputs we're interested in).
	y = scenario.inputs.num_inputs(subinputs, :);

	% Linear interpolation time! (we use extrap so we don't get NaNs back due to rounding errors).
	input = interp1(x, y, t, 'linear', 'extrap');
end
