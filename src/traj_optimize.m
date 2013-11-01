% This function runs a trajectory optimization on the given scenario structure

function [scenario,success] = traj_optimize(scenario, noopt_params)
	disp('Beginning trajectory optimization')
	% Make noopt_params empty, if not given
	if nargin < 2
		noopt_params = [];
	end

	% Just default to the fmincon-based optimizer
	[scenario,success] = traj_opt_fmincon(scenario, noopt_params);

	% Go ahead and generate the numerics (these will probably be needed, and don't take that long to generate).
	scenario = traj_gen_numerical(scenario, noopt_params);
end

% This function generates numerical representations used for analysis of the solution that's been found.
function scenario = traj_gen_numerical(scenario, noopt_params)
	scenario = traj_eval_funcs(scenario, scenario.soln, noopt_params, false);
end
