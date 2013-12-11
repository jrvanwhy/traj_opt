function scenario = traj_eval_funcs(scenario, opt_params, noopt_params)
	disp('Evaluating functions')
	% Go through each phase. Add each function in each phase to the phase as a symbolic variable
	% As an example, this is what generates the states field of phases
	for iter_phase = 1:numel(scenario.phases)
		% Convenience variable
		phase = scenario.phases{iter_phase};

		disp(['	Processing phase ''' phase.names.phase ''''])

		% Parameters for this phase
		phase_params = opt_params(scenario.param_map(1,iter_phase):scenario.param_map(2,iter_phase));

		for iter = 1:numel(phase.functions)
			iter_fcn = phase.functions{iter};
			disp(['		Processing function ''' iter_fcn.name ''''])

			% Call the function using our vecfcn library
			phase.(iter_fcn.name) = opt_eval_vecfcn(iter_fcn.fcn, phase_params);
		end

		% Iterate separately to add time representations, as they require that the duration
		% have already been processed
		for iter = 1:numel(phase.functions)
			iter_fcn = phase.functions{iter};

			% Don't add on the time representation for the duration
			% function
			if strcmp(iter_fcn.name, 'duration')
				continue
			end

			disp(['		Adding time representation for function ''' iter_fcn.name ''''])

			% The name of the corresponding time function
			time_name = ['t_' iter_fcn.name];

			% Add on the time representation. The time spacing depends on the ODE solver type
			switch phase.technique
				case 'midpoint'
					% States are different from everything else
					if strcmp(iter_fcn.name, 'states')
						phase.(time_name) = phase.duration * linspace(0, 1, 1+phase.n_intervals);
					else
						full_grid = linspace(0, 1, 2*phase.n_intervals+1);
						phase.(time_name) = phase.duration * full_grid(2:2:2*phase.n_intervals);
					end

				otherwise
					error(['Technique ''' phase.technique ''' not recognized.'])
			end
		end

		% Copy the phase back over, since phase was produced purely for convenience
		scenario.phases{iter_phase} = phase;
	end
end
