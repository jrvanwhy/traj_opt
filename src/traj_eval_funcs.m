function scenario = traj_eval_funcs(scenario, opt_params, noopt_params)
	disp('Evaluating functions')
	% Go through each phase. Add each function in each phase to the phase as a symbolic variable
	% As an example, this is what generates the states field of phases
	for iter_phase = 1:numel(scenario.phases)
		phase = scenario.phases{iter_phase};

		disp(['	Processing phase ''' phase.names.phase ''''])

		% Parameters for this phase
		phase_params = opt_params(scenario.param_map(1,iter_phase):scenario.param_map(2,iter_phase));

		for iter = 1:numel(phase.functions)
			iter_fcn = phase.functions{iter};
			disp(['		Processing function ''' iter_fcn.name ''''])

			phase.(iter_fcn.name) = iter_fcn.fcn(phase_params);
		end

		scenario.phases{iter_phase} = phase;
	end
end
