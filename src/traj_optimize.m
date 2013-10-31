% This function runs a trajectory optimization on the given scenario structure

function [scenario,success] = traj_optimize(scenario, noopt_params)
	disp('Beginning trajectory optimization')
	% Make noopt_params empty, if not given
	if nargin < 2
		noopt_params = [];
	end

	% Just default to the fmincon-based optimizer
	[scenario,success] = traj_opt_fmincon(scenario, noopt_params);
end
