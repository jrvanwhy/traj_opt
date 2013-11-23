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
			% Error checking
			if ~isempty(scenario_fcn)
				error('Multiple scenario functions given! Only one may be specified')
			end

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

	% Do all the heavy processing on the phases
	scenario = process_phases(scenario);

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

% Handles everything necessary for the scenario function
function scenario = process_scenario_fcn(scenario, scenario_fcn)
	% Check if a scenario function was given. If not, return scenario.
	if isempty(scenario_fcn)
		scenario = scenario;
		return
	end

	disp('	Setting up to call scenario function')

	% Set up the symbolic variables representing the optimization and additional parameters
	opt_params   = sym('x', [scenario.param_maps.n_opt_params   1]);

	% Add the symbolic outputs for each function to the scenario
	scenario = traj_eval_funcs(scenario, opt_params, noopt_params, true);

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
