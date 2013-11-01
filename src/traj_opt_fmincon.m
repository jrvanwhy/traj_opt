% This function runs a numerical optimization using MATLAB's fmincon function

function [scenario,success] = traj_opt_fmincon(scenario, noopt_params)
	% Do the setup, if necessary
	scenario = setup_opt_fmincon(scenario);

	% Run the optimization
	[scenario, success] = run_opt_fmincon(scenario, noopt_params);
end

% This function does all the setup for an fmincon optimization run (if necessary)
function scenario = setup_opt_fmincon(scenario)
	% Exit if this has already been done
	if isfield(scenario, 'fmincon')
		scenario = scenario;
		return
	end

	% Diagnostic message for the user
	disp('	Doing setup for fmincon-based optimization.')

	% Create symbolic variables representing the optimization and non-optimized parameters
	opt_params   = sym('x', [scenario.param_maps.n_opt_params   1]);
	noopt_params = sym('n', [scenario.param_maps.n_noopt_params 1]);

	% Iterate through the cost functions to generate the objective
	objective = 0;
	for iter = 1:numel(scenario.optfcns.costs)
		objective = objective + scenario.optfcns.costs{iter}.fcn(opt_params, noopt_params);
	end

	% Iterate through the constraint functions to generate the constraint
	c   = sym([]);
	ceq = sym([]);
	for iter = 1:numel(scenario.optfcns.constraints)
		% Switch based on the type of constraint
		if scenario.optfcns.constraints{iter}.con_type == '<='
			c   = [c;   reshape(scenario.optfcns.constraints{iter}.fcn(opt_params, noopt_params), [], 1)];
		else
			ceq = [ceq; reshape(scenario.optfcns.constraints{iter}.fcn(opt_params, noopt_params), [], 1)];
		end
	end

	% Clean up symbolic variables
	disp('	Cleaning up symbolic variables')
	syms x n clear
end
