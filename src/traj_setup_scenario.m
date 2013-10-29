% This function does most of the setup and processing for a given scenario.
% It returns a scenario ready for optimization (or generation
% of a fast optimization function, if desired).

% Scenario structure setup:
% interval    A large interval structure containing all the phases.
% name        This scenario's name. Blank if not specified by the user.
% phases      An array of the phase structures in this scenario
% version     The version of the optimizer that generated this scenario
% struct_type The type of this structure -- always 'scenario'

function scenario = traj_setup_scenario(varargin)
	disp('Initializing scenario structure')

	% Basic structure initialization
	scenario.version     = traj_version();
	scenario.struct_type = 'scenario';

	% Initialize the name to a blank string -- it will be changed later
	% if the user specified a name
	scenario.name = '';

	% Also, initialize the scenario function to empty -- it will be filled in if provided
	scenario_fcn = [];

	% Process each input in order
	for iter = 1:numel(varargin)
		% Switch based on the type of this input
		if ischar(varargin{iter})
			% This is a name for the scenario
			disp(['	Processing name ''' varargin{iter} ''''])

			% Check if a name's already been supplied -- if so, error
			if ~isempty(scenario.name)
				error(['Multiple names specified. The first: ''' scenario.name ''', the second: ''' varargin{iter} ''''])
			end

			% Add the name to the scenario
			scenario.name = varargin{iter};
		elseif isa(varargin{iter}, 'function_handle')
			% This is the scenario function
			% Save it to a variable for later use.
			scenario_fcn = varargin{iter};

			% Verify that the number of output arguments is positive
			if (nargout(scenario_fcn) < 1)
				error('The scenario function must have a constant positive number of output arguments')
			end

		elseif isstruct(varargin{iter}) && isfield(varargin{iter}, 'struct_type') && strcmp(varargin{iter}.struct_type, 'dynamic phase')
			% This is a new phase.
			disp(['	Processing phase ''' varargin{iter}.names.phase ''''])

			% Update the version of this input struct
			varargin{iter} = traj_version_update(varargin{iter});

			% Call the phase processing function
			scenario = add_phase(scenario, varargin{iter});
		else
			% Spit out an error -- this input is not recognized
			error(['Input ' num2str(iter) ' not recognized'])
		end
	end

	% Make sure at least one phase was specified.
	if ~isfield(scenario, 'phases')
		error('At least one phase must be passed!')
	end

	% Map the optimization and non-optimized parameters to the interval's
	% various parameter sets
	scenario.param_maps = map_interval_params(scenario);

	% Call and process data from the scenario function.
	scenario = process_scenario_fcn(scenario, scenario_fcn);

	% Generate the final objecttive and constraint functions
	scenario.optfcns = gen_optfcns(scenario);
end

% This function adds a phase to the given scenario
function scenario = add_phase(scenario, phase)
	% Append to the phases list (creating it if necessary)
	if isfield(scenario, 'phases')
		scenario.phases(end+1) = phase;
	else
		scenario.phases(1) = phase;
	end
end

% This function generates the optfcns structure, containing the optimization costs and constraints
function optfcns = gen_optfcns(scenario)
	% Initialize as empty, so we can simply append to them.
	optfcns.constraints = [];
	optfcns.costs       = [];

	% Append scenario function costs and constraints
	if isfield(scenario.scenfun, 'costs')
		optfcns.costs = [optfcns.costs, scenario.scenfun.costs];
	end
	if isfield(scenario.scenfun, 'constraints')
		optfcns.constraints = [optfcns.constraints, scenario.scenfun.constraints];
	end

	% Append phase costs and constraints
	if isfield(scenario.phases.phase_interval, 'costs')
		optfcns.costs = [optfcns.costs, scenario.phases.phase_interval.costs];
	end
	if isfield(scenario.phases.phase_interval, 'constraints')
		optfcns.constraints = [optfcns.constraints, scenario.phases.phase_interval.constraints];
	end
end

% This function creates functions mapping the optimization parameters to the
% interval's parameters
function param_maps = map_interval_params(scenario)
	% Count the optimization and non-optimization parameters
	param_maps.n_opt_params = numel(scenario.phases.phase_interval.start_params)  + ...
	                          numel(scenario.phases.phase_interval.end_params)    + ...
	                          numel(scenario.phases.phase_interval.int_params)    + ...
	                          numel(scenario.phases.phase_interval.shared_params) + ...
	                          numel(scenario.phases);
	param_maps.n_noopt_params = numel(scenario.phases.phase_interval.noopt_params);

	% Set up symbolic variables representing the optimization and additional parameters
	disp('	Creating optimization and non-optimized symbolic variables')
	opt_params   = sym('x', [param_maps.n_opt_params,   1]);
	noopt_params = sym('n', [param_maps.n_noopt_params, 1]);

	% Positions of the interval's param sets within opt_params
	start_params_pos   = [0 0];
	end_params_pos     = [0 0];
	int_params_pos     = [0 0];
	shared_params_pos  = [0 0];
	duration_param_pos = 0;
	noopt_params_pos   = [0 0];

	% Go through each phase, adding parameter map functions
	for iter_phase = 1:numel(scenario.phases)
		% Update positions in parameter list
		% We stagger each of them in order per phase.
		start_params_pos(1)  = duration_param_pos   + 1;
		start_params_pos(2)  = start_params_pos(1)  + numel(scenario.phases(iter_phase).phase_interval.start_params)  - 1;
		end_params_pos(1)    = start_params_pos(2)  + 1;
		end_params_pos(2)    = end_params_pos(1)    + numel(scenario.phases(iter_phase).phase_interval.end_params)    - 1;
		int_params_pos(1)    = end_params_pos(2)    + 1;
		int_params_pos(2)    = int_params_pos(1)    + numel(scenario.phases(iter_phase).phase_interval.int_params)    - 1;
		shared_params_pos(1) = int_params_pos(2)    + 1;
		shared_params_pos(2) = shared_params_pos(1) + numel(scenario.phases(iter_phase).phase_interval.shared_params) - 1;
		duration_param_pos   = shared_params_pos(2) + 1;

		% Noopt parameters are special, because they're a different argument
		noopt_params_pos(1)  = noopt_params_pos(2) + 1;
		noopt_params_pos(2)  = noopt_params_pos(1) + numel(scenario.phases(iter_phase).phase_interval.noopt_params)  - 1;

		% Generate the functions
		disp(['		Generating param map functions for phase ''' scenario.phases(iter_phase).names.phase '''']);
		param_maps.start_params{iter_phase}  = matlabFunction(...
			opt_params(start_params_pos(1) :start_params_pos(2) ), 'vars', {opt_params, noopt_params});
		param_maps.end_params{iter_phase}    = matlabFunction(...
			opt_params(end_params_pos(1)   :end_params_pos(2)   ), 'vars', {opt_params, noopt_params});
		param_maps.int_params{iter_phase}    = matlabFunction(...
			opt_params(int_params_pos(1)   :int_params_pos(2)   ), 'vars', {opt_params, noopt_params});
		param_maps.shared_params{iter_phase} = matlabFunction(...
			opt_params(shared_params_pos(1):shared_params_pos(2)), 'vars', {opt_params, noopt_params});
		param_maps.noopt_params{iter_phase}  = matlabFunction(...
			noopt_params,                                          'vars', {opt_params, noopt_params});
		param_maps.duration_param{iter_phase} = matlabFunction(...
			opt_params(duration_param_pos),                        'vars', {opt_params, noopt_params});
	end

	% TODO: Correctly deal with noopt_params, making use of the parameter names.

	% Clean up
	disp('	Cleaning up symbolic variables')
	syms x n clear
end

% Handles everything necessary for the scenario function
function scenario = process_scenario_fcn(scenario, scenario_fcn)
	% Check if a scenario function was given. If not, return scenario and exit.
	if isempty(scenario_fcn)
		scenario = scenario;
		return
	end

	disp('	Setting up to call scenario function')

	% Set up the symbolic variables representing the optimization and additional parameters
	opt_params   = sym('x', [scenario.param_maps.n_opt_params   1]);
	noopt_params = sym('n', [scenario.param_maps.n_noopt_params 1]);

	% Go through each phase. Add each function in each phase to the phase as a symbolic variable
	% As an example, this is what generates the states field of phases
	for iter_phase = 1:numel(scenario.phases)
		fcn_names = fieldnames(scenario.phases(iter_phase).phase_interval.funcs);
		for iter = 1:numel(fcn_names)
			iter_fcn = fcn_names{iter};
			disp(['		Processing function ''' iter_fcn ''''])

			% Add in the field, so it may be appended to. Note that this must be set to the right type
			% at initialization, so we'll make it symbolic from the start.
			scenario.phases(iter_phase).(iter_fcn) = sym([]);

			% Iterate through all sub-functions in the array, appending them to the new field
			for iter_col = 1:numel(scenario.phases(iter_phase).phase_interval.funcs.(iter_fcn))
				% Grab the function for this column
				col_fcn = scenario.phases(iter_phase).phase_interval.funcs.(iter_fcn){iter_col};

				% Call the column function
				scenario.phases(iter_phase).(iter_fcn)(:,end+1) = simplify(col_fcn(...
					scenario.param_maps.start_params{iter_phase}(opt_params, noopt_params),  ...
					scenario.param_maps.end_params{iter_phase}(opt_params, noopt_params),    ...
					scenario.param_maps.int_params{iter_phase}(opt_params, noopt_params),    ...
					scenario.param_maps.shared_params{iter_phase}(opt_params, noopt_params), ...
					scenario.param_maps.noopt_params{iter_phase}(opt_params, noopt_params),  ...
					scenario.param_maps.duration_param{iter_phase}(opt_params, noopt_params)));
			end
		end
	end

	% Call the scenario function, capturing all output arguments
	disp('	Calling scenario function')
	[scenario_varouts{1:nargout(scenario_fcn)}] = scenario_fcn(scenario);

	% Initialize scenario function fields
	scenario.scenfun.constraints = {};
	scenario.scenfun.costs       = {};

	% Iterate through the costs and constraints, adding them to phases
	disp('	Adding costs and constraints to scenario function representation.')
	for idx = 1:numel(scenario_varouts)
		% This is the structure returned in this output
		output = scenario_varouts{idx};

		% Check that this output has a valid type string
		if ~isfield(output, 'struct_type')
			% Spit out an error message; this data type is invalid
			error(['Scenario system function output ' num2str(idx+1) ' is not of a valid data type'])
		end

		% Use a switch statement to decide how to handle this output
		switch output.struct_type
			case 'constraint'
				% Let the user know things are happening
				disp(['		Processing constraint ''' output.name ''''])

				% Convert the constraint function to an anonymous function, then append it to the
				% constraints field.
				output.fcn = matlabFunction(simplify(output.fcn), 'vars', ...
					{opt_params, noopt_params});
				scenario.scenfun.constraints{end+1} = output;

			case 'cost'
				% Progress message so the user doesn't think it's frozen
				disp(['		Processing cost ''' output.name ''''])

				% Convert the cost function to an anonymous function with appropriate inputs, then
				% append it to the costs field
				output.fcn = matlabFunction(output.fcn, 'vars', ...
					{opt_params, noopt_params});
				scenario.scenfun.costs{end+1} = output;

			otherwise
				% Spit out an error -- this is not an acceptable structure type
				error(['Structure type ''' output.struct_type ''' not a valid return type for the dynamic system function.'])
		end
	end

	disp(['	Cleaning up symbolic variables'])
end