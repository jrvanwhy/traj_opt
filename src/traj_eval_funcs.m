function scenario = traj_eval_funcs(scenario, opt_params, noopt_params, symbolic)
	% Go through each phase. Add each function in each phase to the phase as a symbolic variable
	% As an example, this is what generates the states field of phases
	for iter_phase = 1:numel(scenario.phases)
		fcn_names = fieldnames(scenario.phases(iter_phase).phase_interval.funcs);
		for iter = 1:numel(fcn_names)
			iter_fcn = fcn_names{iter};
			disp(['		Processing function ''' iter_fcn ''''])

			% Add in the field, so it may be appended to. Note that this must be set to the right type
			% at initialization, so we'll make it symbolic from the start.
			if symbolic
				scenario.phases(iter_phase).(iter_fcn) = sym([]);
			else
				scenario.phases(iter_phase).(iter_fcn) = [];
			end

			% Iterate through all sub-functions in the array, appending them to the new field
			for iter_col = 1:numel(scenario.phases(iter_phase).phase_interval.funcs.(iter_fcn))
				% Grab the function for this column
				col_fcn = scenario.phases(iter_phase).phase_interval.funcs.(iter_fcn){iter_col};

				% Call the column function
				scenario.phases(iter_phase).(iter_fcn)(:,end+1) = col_fcn(...
					scenario.param_maps.start_params{iter_phase}(opt_params, noopt_params),  ...
					scenario.param_maps.end_params{iter_phase}(opt_params, noopt_params),    ...
					scenario.param_maps.int_params{iter_phase}(opt_params, noopt_params),    ...
					scenario.param_maps.shared_params{iter_phase}(opt_params, noopt_params), ...
					scenario.param_maps.noopt_params{iter_phase}(opt_params, noopt_params),  ...
					scenario.param_maps.duration_param{iter_phase}(opt_params, noopt_params));
			end
		end
	end
end
