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

	% Call and process data from the scenario function.
	scenario = process_scenario_fcn(scenario, scenario_fcn);
end

% This function adds a phase to the given scenario
function scenario = add_phase(scenario, phase)
	% Append to the phases list (create it if necessary)
	% Also, append the intervals (or just copy in interval if none exist yet)
	if isfield(scenario, 'phases')
		scenario.phases(end+1) = phase;
		scenario.interval = concatenate_intervals(scenario.interval, phase.interval);
	else
		scenario.phases(1) = phase;
		scenario.interval  = phase.interval;
	end
end

% This function generates the states and phases for the given scenario
function scenario = gen_states_inputs(scenario)
	error('TODO: this')
end

% This function concatenates two intervals, producing one larger interval
function intout = concatenate_intervals(intin1, intin2)
	error('TODO: this')
end

% Handles everything necessary for the scenario function
function scenario = process_scenario_fcn(scenario, scenario_fcn)
	% Check if a scenario function was given. If not, return scenario and exit.
	if isempty(scenario_fcn)
		scenario = scenario;
		return
	end

	% Set up the symbolic variables representing the optimization and additional parameters
	opt_params   = sym('x', [numel(scenario.phases.phase_interval.start_params)  + ...
	                         numel(scenario.phases.phase_interval.end_params)    + ...
	                         numel(scenario.phases.phase_interval.int_params)    + ...
	                         numel(scenario.phases.phase_interval.shared_params) + 1, 1]); % Add up all parameters (including 1 extra for duration).
	noopt_params = sym('n', [numel(scenario.phases.phase_interval.noopt_params) 1]);

	% Variables for keeping track of our position in the parameter list
	start_params_pos   = [0 0];
	end_params_pos     = [0 0];
	int_params_pos     = [0 0];
	shared_params_pos  = [0 0];
	noopt_params_pos   = [0 0];
	duration_param_pos = 0;

	% Go through each phase. Add each function in each phase to the phase as a symbolic variable
	% As an example, this is what generates the states field of phases
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
					opt_params(start_params_pos(1):start_params_pos(2)), ...
					opt_params(end_params_pos(1)  : end_params_pos(2)),  ...
					opt_params(int_params_pos(1)  : int_params_pos(2)),  ...
					opt_params(shared_params_pos(1) : shared_params_pos(2)), ...
					noopt_params(noopt_params_pos(1) : noopt_params_pos(2)), ...
					opt_params(duration_param_pos)));
			end
		end
	end

	% Call the scenario function, capturing all output arguments
	[scenario_varouts{1:nargout(scenario_fcn)}] = scenario_fcn(scenario);
end
