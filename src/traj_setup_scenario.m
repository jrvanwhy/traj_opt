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

		elseif isstruct(varargin{iter}) && isfield(varargin{iter}, 'struct_type') && strcmp(varargin{iter}.struct_type, 'dynamic phase')
			% This is a new phase.
			disp(['	Processing phase ''' varargin{iter}.names.phase ''''])

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
	error('TODO: this')
end
