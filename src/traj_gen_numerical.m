% This generates numerical forms for many of the values in the scenario.
% Manipulating numerical values after an optimization is done is faster and often more convenient
% than manipulating symbolic values.

function scenario = traj_gen_numerical(scenario)
	sym_params                  = [scenario.params.params; scenario.add_params.params];
	sym_param_vals              = [scenario.minimum(:,end); scenario.add_params.value];
	scenario.num_duration       = double(subs(scenario.duration,       sym_params, sym_param_vals));
	scenario.inputs.num_inputs  = double(subs(scenario.inputs.inputs,  sym_params, sym_param_vals));
	scenario.states.num_dstates = double(subs(scenario.states.dstates, sym_params, sym_param_vals));
	scenario.states.num_states  = double(subs(scenario.states.states,  sym_params, sym_param_vals));
end
