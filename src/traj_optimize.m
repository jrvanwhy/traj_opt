% This function runs a trajectory optimization on the given scenario structure

function [scenario,success] = traj_optimize(scenario, noopt_params, iter_fcn)
	disp('Beginning trajectory optimization')
	% Make noopt_params empty, if not given
	if nargin < 2
		noopt_params = [];
	end

	% Make iter_fcn empty, if not given
	if nargin < 3
		iter_fcn = [];
	end

	% Iteration function to pass to the optimization function
	if isempty(iter_fcn)
		opt_iter_fcn = [];
	else
		opt_iter_fcn = @opt_iterfcn;
	end

	% Just default to the fmincon-based optimizer
	[scenario,success] = traj_opt_fmincon(scenario, noopt_params, opt_iter_fcn);

	% Go ahead and generate the numerics (these will probably be needed, and don't take that long to generate).
	scenario = traj_gen_numerical(scenario, noopt_params);

	% Optimizer iteration function; this evaluates numerical values
	% then calls the user's iteration function
	function stop = opt_iterfcn(params, optimValues, state)
		scenario = traj_eval_funcs(scenario, params, noopt_params);
		iter_fcn(scenario);

		% Set our output to continue the optimization
		stop = false;
	end
end

% This function generates numerical representations used for analysis of the solution that's been found.
function scenario = traj_gen_numerical(scenario, noopt_params)
	disp('	Generating numerical values')
	scenario = traj_eval_funcs(scenario, scenario.soln, []);
end
