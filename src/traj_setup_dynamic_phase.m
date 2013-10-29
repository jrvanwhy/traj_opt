% This function creates a new dynamic phase structure.
% Parameters:
%     name         The name for this phase (optional -- leave empty for an empty string).
%     dynsys_fcn   A function of the form [dx,...] = dynsys_fcn(x, u, add_params, noopt_params)
%                  describing the dynamics, costs, and constraints for this phase. If state, input, add_params, or noopt_params
%                  are empty, those should be omitted from the parameter list for dynsys_fcn.
%     state        A cell array of strings naming each substate of the state for this phase.
%     input        A cell array of strings naming each subinput of the input for this phase.
%     n_intervals  The number of time intervals in the discrete approximation to this phase.
%                  This is roughly equivalent to the number of timesteps for an ODE solver.
%     add_params   Optimization parameters specific to this phase that are not the state or the input (cell array of names).
%     noopt_params Parameters for this phase that should not be affected by the optimizer, but will be set
%                  by the user prior to running the optimization (cell array of names).

% The phases structure has the following fields:
%     duration    The duration, for the scenario function and later processing
%     dynsys      A structure containing anonymous functions derived from the dynsys_fcn argument
%         constraints A cell array of constraints
%         costs       A cell array of cost functions
%         dx          The dynamics function (ẋ = f(x, u, a, n))
%     inputs      This contains the inputs for this phase (for the scenario function and later processing).
%     interval    An interval structure representing one interval within the dynamic phase. Each function inside is defined
%                 in terms of the start params, end params, int_params, shared params, noopt_params, and duration (in that order).
%         name          This interval type's name
%         costs         The per-interval costs for this interval type
%         constraints   A cell array of constraints for this interval type (per-interval)
%         start_funcs   Additional functions to be evaluated only at the start of the phase
%         funcs         Miscellaneous functions (state, duration, input, ...). Most map parameters to more useful values.
%         start_params  Parameters shared with the "previous" interval of this type.
%         end_params    Parameters shared with the "next" (dynamically) interval of this type.
%         shared_params Parameters shared by all intervals of this type
%         noopt_params  Shared parameters not optimized by the optimizer
%         int_params    Parameters specific to each interval of this type
%     n_intervals The number of intervals in this phase
%     names       A structure containing names for this phase and various sub-fields. It is as follows:
%         add_params   Cell array of additional parameter names
%         input        Cell array of input names
%         noopt_params Cell array of names for the non-optimized parameters.
%         phase        The name of this phase
%         state        The cell array of state names
%     phase_interval An interval structure representing the entire phase (cloned from the member 'interval' above)
%     states      The states for this phase (for the scenario function and later processing).
%     technique   A string indicating the trajectory optimization technique used, such as 'dircol 1'
%     struct_type      The type of this struct -- always 'dynamic phase'
%     version    The version of the optimizer this phase was last generated by

function phase = traj_setup_dynamic_phase(name, dynsys_fcn, state, input, n_intervals, add_params, noopt_params)
	% This section of this function sets default values for its arguments
	% The default value for name is the empty string
	if nargin < 1 || isempty(name)
		name = '';
	end

	% The default value for state and input are empty cell arrays.
	if nargin < 3 || isempty(state)
		state = {};
	end
	if nargin < 4 || isempty(input)
		input = {};
	end

	% The default value for n_intervals is 30
	if nargin < 5 || isempty(n_intervals)
		n_intervals = 30;
	end

	% The default value for add_params and noopt_params are empty cell arrays
	if nargin < 6 || isempty(add_params)
		add_params = {};
	end
	if nargin < 7 || isempty(noopt_params)
		noopt_params = {};
	end

	% This is the input sanity verification section
	% name should be a character array (string)
	if ~ischar(name)
	;	error('name must be empty or a string')
	end

	% dynsys_fcn must be a function handle with a fixed number of arguments and return values
	if ~isa(dynsys_fcn, 'function_handle') || nargout(dynsys_fcn) < 1
		error('dynsys_fcn must be a function handle with at least one output arguments (and no varargout)')
	end
	% Check if dynsys_fcn has the right number of input arguments
	if nargin(dynsys_fcn) ~= ~isempty(state) + ~isempty(input) + ~isempty(add_params) + ~isempty(noopt_params)
		error('dynsys_fcn does not have the correct number of arguments.')
	end

	% State and input must be cell arrays of strings
	if ~iscellstr(state)
		error('state must be a cell array of strings')
	end
	if ~iscellstr(input)
		error('input must be a cell array of strings')
	end

	% Make sure the number of intervals makes sense
	% The mod(n_intervals, 1) boolean value is equivalent to mod(n_intervals, 1) ~= 0, which is only true
	% if n_intervals is not an integer.
	if mod(n_intervals, 1) || n_intervals < 0
		error('n_intervals must be a nonnegative integer')
	end

	% add_params and noopt_params must be cell arrays of strings
	if ~iscellstr(add_params)
		error('add_params must be a cell array of strings')
	end
	if ~iscellstr(noopt_params)
		error('noopt_params must be a cell array of strings')
	end

	% Default parameters and input sanity checks done. The actual core of the function is below.

	% Spit out a simple message to indicate stuff's happening
	disp(['Beginning setup for phase ''' name ''''])

	% Set struct version and type
	phase.version     = traj_version();
	phase.struct_type = 'dynamic phase';

	% Copy the name and number of intervals into the phase structure
	phase.n_intervals = n_intervals;
	phase.names.phase = name;

	% Copy over state, input, and parameter name arrays (converting to column vectors where necessary)
	phase.names.add_params   = add_params(:);
	phase.names.input        = input(:);
	phase.names.noopt_params = noopt_params(:);
	phase.names.state        = state(:);

	% Generate dynamic system function-related functions in phase
	phase = gen_dynsys_fcns(phase, dynsys_fcn);

	% Generate the interval structure
	phase = gen_interval(phase);

	% Generate the phase's interval structure
	phase.phase_interval = clone_interval(phase.interval, phase.n_intervals);

	% Add on the duration function (because it's special and important)
	%phase = phase_gen_duration(phase); % TODO: This

	% Let the user know when this function terminates
	disp(['Setup for phase ''' name ''' completed.'])
end

% This function 'clones' an interval, generating a larger interval that consisting of multiple copies
% of the input interval combined into one.
function int_out = clone_interval(int_in, count)
	% Informative message for the user
	disp(['	Cloning interval ''' int_in.name ''' ' num2str(count) ' times.'])

	% Copy name
	int_out.name = int_in.name;

	% Define starting and ending parameters (they are simply numbered versions of the input interval's
	% starting and ending parameters)
	for iter = 1:numel(int_in.start_params)
		int_out.start_params{iter,1} = [int_in.start_params{iter} '/1'];
	end
	for iter = 1:numel(int_in.end_params)
		int_out.end_params{iter,1} = [int_in.end_params{iter} '/' num2str(count+1)];
	end

	% The shared parameters are simply copied over, as are the noopt parameters
	int_out.shared_params = int_in.shared_params(:);
	int_out.noopt_params  = int_in.noopt_params(:);

	% The internal states become interior params; the inputs stay interior params but are duplicated.
	% We'll insert the inputs first, then the states
	int_out.int_params = [];
	for iter1 = 1:count
		% Append numbers to the inputs
		for iter2 = 1:numel(int_in.int_params)
			int_out.int_params{end+1,1} = [int_in.int_params{iter2} '/' num2str(iter1)];
		end
	end
	for iter1 = 2:count
		% Append numbers to each sub-state
		for iter2 = 1:numel(int_in.end_params)
			int_out.int_params{end+1,1} = [int_in.end_params{iter2} '/' num2str(iter1)];
		end
	end

	% Generate symbolic variables for each input for the new interval's functions
	sym_start_params  = sym('b', size(int_out.start_params));
	sym_end_params    = sym('e', size(int_out.end_params));
	sym_int_params    = sym('i', size(int_out.int_params));
	sym_shared_params = sym('s', size(int_out.shared_params));
	sym_noopt_params  = sym('n', size(int_out.noopt_params));
	sym_duration      = sym('t', 'real');

	% Create sub-interval input representations
	disp('		Creating sub-interval input representations')
	subint_start_params = [sym_start_params, ...
		reshape(sym_int_params(count*numel(int_in.int_params)+1:...
			count*(numel(int_in.int_params)+numel(int_in.end_params))-numel(int_in.end_params)), ...
			numel(int_in.end_params), count-1)];
	subint_end_params = [subint_start_params(:,2:end), sym_end_params];
	subint_int_params = reshape(sym_int_params(1:count*numel(int_in.int_params)), numel(int_in.int_params), count);

	% Clone the constraints if there are any
	if isfield(int_in, 'constraints')
		disp('		Cloning constraints')
		int_out.constraints = [];
		for iter1 = 1:numel(int_in.constraints)
			disp(['			Processing constraint ''' int_in.constraints{iter1}.name ''''])
			for iter2 = 1:count
				% Copy over the basic constraint structure
				int_out.constraints{end+1} = int_in.constraints{iter1};

				% Append a number to the constraint's name
				int_out.constraints{end}.name = [int_in.constraints{iter1}.name '/' num2str(iter2)];

				% Update the constraint function to reflect the new function inputs
				int_out.constraints{end}.fcn = int_in.constraints{iter1}.fcn(...
					subint_start_params(:,iter2), ...
					subint_end_params(:,iter2),   ...
					subint_int_params(:,iter2),   ...
					sym_shared_params,            ...
					sym_noopt_params,             ...
					sym_duration/count);

				% Convert the function to a symbolic function
				int_out.constraints{end}.fcn = matlabFunction(int_out.constraints{end}.fcn, 'vars', ...
					{sym_start_params, sym_end_params, sym_int_params, sym_shared_params, sym_noopt_params, sym_duration});
			end
		end
	end

	% Clone the costs, if any exist
	if isfield(int_in, 'costs')
		disp('		Cloning costs')
		int_out.costs = [];
		for iter1 = 1:numel(int_in.costs)
			disp(['			Processing cost ''' int_in.costs{iter1}.name ''''])
			for iter2 = 1:count
				% Copy over the basic cost structure
				int_out.costs{end+1} = int_in.costs{iter1};

				% Append a number to the cost's name
				int_out.costs{end}.name = [int_in.costs{iter1}.name '/' num2str(iter2)];

				% Update the cost function to reflect the new function inputs
				int_out.costs{end}.fcn = int_in.costs{iter1}.fcn(...
					subint_start_params(:,iter2), ...
					subint_end_params(:,iter2),   ...
					subint_int_params(:,iter2),   ...
					sym_shared_params,            ...
					sym_noopt_params,             ...
					sym_duration/count);

				% Convert the function to a symbolic function
				int_out.costs{end}.fcn = matlabFunction(int_out.costs{end}.fcn, 'vars', ...
					{sym_start_params, sym_end_params, sym_int_params, sym_shared_params, sym_noopt_params, sym_duration});
			end
		end
	end

	% Iterate through starting functions, generating them with the correct parameters
	if isfield(int_in, 'start_funcs')
		% The list of function names
		fcn_names = fieldnames(int_in.start_funcs);

		% Iterate through them
		for iter = 1:numel(fcn_names)
			disp(['		Processing starting function ''' fcn_names{iter} ''''])

			% Call it with the correct parameters first and simplify
			int_out.funcs.(fcn_names{iter}){1,1} = simplify(int_in.start_funcs.(fcn_names{iter})(...
				subint_start_params(:,1), ...
				subint_end_params(:,1),   ...
				subint_int_params(:,1),   ...
				sym_shared_params,        ...
				sym_noopt_params,         ...
				sym_duration/count));
		end
	end

	% Iterate through the additional functions, generating new ones that that return their value for each sub-interval
	if isfield(int_in, 'funcs')
		% The list of function names
		fcn_names = fieldnames(int_in.funcs);

		for iter = 1:numel(fcn_names)
			disp(['		Processing function ''' fcn_names{iter} ''''])

			% If the output function is empty, initialize it
			if (~isfield(int_out.funcs, fcn_names{iter}))
				int_out.funcs.(fcn_names{iter}) = {};
			end

			% Iterate through each sub-interval, creating a new function for each sub-interval
			for iter2 = 1:count
				% Create the symbolic form of the new function
				int_out.funcs.(fcn_names{iter}){1,end+1} = simplify(int_in.funcs.(fcn_names{iter})(...
					subint_start_params(:,iter2), ...
					subint_end_params(:,iter2),   ...
					subint_int_params(:,iter2),   ...
					sym_shared_params,            ...
					sym_noopt_params,             ...
					sym_duration/count));
			end
		end
	end

	% Convert additional functions to anonymous functions
	if isfield(int_out, 'funcs')
		% Iterate through each function
		fcn_names = fieldnames(int_out.funcs);

		for iter = 1:numel(fcn_names)
			disp(['		Converting function ''' fcn_names{iter} ''''])

			% Iterate through each element, converting one-by-one for performance.
			for iter2 = 1:numel(int_out.funcs.(fcn_names{iter}))
				% Do the conversion
				int_out.funcs.(fcn_names{iter}){iter2} = matlabFunction(...
					int_out.funcs.(fcn_names{iter}){iter2}, 'vars', ...
					{sym_start_params, sym_end_params, sym_int_params, sym_shared_params, sym_noopt_params, sym_duration});
			end
		end
	end

	% Clean up the symbolic variables
	disp('		Cleaning up symbolic variables')
	syms b e i s n t clear
end

% This function generates anonymous functions representing the dynamic system function
% and places them into phase
function phase = gen_dynsys_fcns(phase, dynsys_fcn)
	% Diagnostic message to inform the user of the framework's progress
	disp('	Generating dynamic system function representations')

	% Create symbolic variables representing the various parameters to the dynamic system function
	% We use these symbolic variables to create representations of the dynamic system functions as anonymous function.
	sym_state        = sym('x', [numel(phase.names.state)        1]);
	sym_input        = sym('u', [numel(phase.names.input)        1]);
	sym_add_params   = sym('a', [numel(phase.names.add_params)   1]);
	sym_noopt_params = sym('n', [numel(phase.names.noopt_params) 1]);

	% Create a cell array with the parameters to the dynamic system function. This allows us to call it with the
	% correct number of arguments.
	dynsys_args = {};
	if ~isempty(phase.names.state)
		dynsys_args = [dynsys_args, {sym_state}];
	end
	if ~isempty(phase.names.input)
		dynsys_args = [dynsys_args, {sym_input}];
	end
	if ~isempty(phase.names.add_params)
		dynsys_args = [dynsys_args, {sym_add_params}];
	end
	if ~isempty(phase.names.noopt_params)
		dynsys_args = [dynsys_args, {sym_noopt_params}];
	end

	% Call the dynamic system function, retrieving all outputs
	[dynsys_dx,dynsys_varouts{1,1:nargout(dynsys_fcn)-1}] = dynsys_fcn(dynsys_args{:});

	% Set up phases.dynsys.dx, simplifying it in the process for performance
	disp('	Generating dynamics ODE function')
	phase.dynsys.dx = matlabFunction(simplify(dynsys_dx), 'vars', {sym_state, sym_input, sym_add_params, sym_noopt_params});

	% Initialize dynamic system fields
	phase.dynsys.constraints = {};
	phase.dynsys.costs       = {};

	% Iterate through the costs and constraints, adding them to phases
	disp('	Adding costs and constraints to dynamic system representation.')
	for idx = 1:numel(dynsys_varouts)
		% This is the structure returned in this output
		output = dynsys_varouts{idx};

		% Check that this output has a valid type string
		if ~isfield(output, 'struct_type')
			% Spit out an error message; this data type is invalid
			error(['Dynamic system function output ' num2str(idx+1) ' is not of a valid data type'])
		end

		% Use a switch statement to decide how to handle this output
		switch output.struct_type
			case 'constraint'
				% Let the user know things are happening
				disp(['		Processing constraint ''' output.name ''''])

				% Convert the constraint function to an anonymous function, then append it to the
				% constraints field.
				output.fcn = matlabFunction(simplify(output.fcn), 'vars', ...
					{sym_state, sym_input, sym_add_params, sym_noopt_params});
				phase.dynsys.constraints = [phase.dynsys.constraints {output}];

			case 'cost'
				% Progress message so the user doesn't think it's frozen
				disp(['		Processing cost ''' output.name ''''])

				% Convert the cost function to an anonymous function with appropriate inputs, then
				% append it to the costs field
				output.fcn = matlabFunction(output, 'vars', ...
					{sym_state, sym_input, sym_add_params, sym_noopt_params});
				phase.dynsys.costs{end+1} = output;

			otherwise
				% Spit out an error -- this is not an acceptable structure type
				error(['Structure type ''' output.struct_type ''' not a valid return type for the dynamic system function.'])
		end
	end

	% Clean up our symbolic variables
	syms x u a n clear
end

% This generates the interval structure for each step in the dynamic simulation of this phase
function phase = gen_interval(phase)
	% Currently, only one trajectory optimization technique ('dircol 1') is supported.
	phase.technique = 'dircol 1';

	% Go ahead and call the specific technique's interval generation function.
	phase = gen_int_dircol_1(phase);
end

% This generates the interval for a basic first-order direct collocation optimization.
% This treats the state as piecewise linear, linearly interpolating it between the state inputs.
% It then samples the state at the midpoint of the interval, evaluates the dynamics at that location,
% then compares that with the state's derivative (according to its slope) and constrains the two derivative estimates
% to match. Simple midpoint integration is used to integrate the cost function.
function phase = gen_int_dircol_1(phase)
	% Informative message for the user
	disp('	Optimizing using 1st-order direct collocation strategy.')

	% The shared parameters are just the add_params specified by the user; similarly,
	% the noopt_params are just the user-specified noopt params.
	phase.interval.shared_params = [phase.names.add_params];
	phase.interval.noopt_params  = phase.names.noopt_params;

	% The interval's name is the phase's name
	phase.interval.name = phase.names.phase;

	% The start and end params are just the state
	phase.interval.start_params = phase.names.state;
	phase.interval.end_params   = phase.names.state;

	% The only internal parameters are the inputs
	phase.interval.int_params = phase.names.input;

	% Create symbolic representations of the various parameters
	sym_start_params  = sym('b', size(phase.interval.start_params));
	sym_end_params    = sym('e', size(phase.interval.end_params));
	sym_int_params    = sym('i', size(phase.interval.int_params));
	sym_shared_params = sym('s', size(phase.interval.shared_params));
	sym_noopt_params  = sym('n', size(phase.interval.noopt_params));
	sym_duration      = sym('t', 'real');

	% Generate the collocation constraint. Note that I've multiplied both sides by the duration of the interval,
	% so as to avoid dividing by the duration in the constraint definition.
	disp('		Creating collocation constraint')
	phase.interval.constraints{1} = traj_create_constraint('collocation', ...
		sym_end_params - sym_start_params, '=', ...
		phase.dynsys.dx((sym_start_params + sym_end_params)/2, sym_int_params, ...
			sym_shared_params, sym_noopt_params) * sym_duration);

	% Convert the user's constraints to the interval's function form
	disp('		Converting user-supplied constraints to the interval''s function form')
	for iter = 1:numel(phase.dynsys.constraints)
		disp(['			Processing constraint ''' phase.dynsys.constraints{iter}.name ''''])
		phase.interval.constraints{iter+1}     = phase.dynsys.constraints{iter};
		phase.interval.constraints{iter+1}.fcn = phase.dynsys.constraints{iter}.fcn(...
			(sym_start_params + sym_end_params)/2, ...
			sym_int_params, ...
			sym_shared_params, ...
			sym_noopt_params);
	end

	% Convert all constraints to function handles
	disp('		Converting interval constraints to anonymous functions')
	for iter = 1:numel(phase.interval.constraints)
		disp(['			Processing constraint ''' phase.interval.constraints{iter}.name ''''])
		phase.interval.constraints{iter}.fcn = matlabFunction(phase.interval.constraints{iter}.fcn, 'vars', ...
			{sym_start_params, sym_end_params, sym_int_params, sym_shared_params, sym_noopt_params, sym_duration});
	end

	% Like the constraints, convert the user's costs to the interval's function form then into function handles.
	disp('		Processing user costs into interval constraints')
	for iter = 1:numel(phase.dynsys.costs)
		disp(['			Processing cost ''' phase.dynsys.costs{iter}.name ''''])
		phase.interval.costs{iter} = phase.dynsys.costs{iter};
		phase.interval.costs{iter}.fcn = phase.dynsys.costs{iter}.fcn(...
			(sym_start_params + sym_end_params)/2, ...
			sym_int_params, ...
			sym_shared_params(2:end), ...
			sym_noopt_params);
		phase.interval.costs{iter}.fcn = matlabFunction(phase.interval.costs{iter}.fcn, 'vars', ...
			{sym_start_params, sym_end_params, sym_int_params, sym_shared_params, sym_noopt_params, sym_duration});
	end

	% Put in useful functions
	phase.interval.funcs.duration = matlabFunction(sym_duration, 'vars', ...
		{sym_start_params, sym_end_params, sym_int_params, sym_shared_params, sym_noopt_params, sym_duration});
	% There is one more state than there are intervals...
	phase.interval.start_funcs.states = matlabFunction(sym_start_params, 'vars', ...
		{sym_start_params, sym_end_params, sym_int_params, sym_shared_params, sym_noopt_params, sym_duration});
	phase.interval.funcs.states       = matlabFunction(sym_end_params, 'vars', ...
		{sym_start_params, sym_end_params, sym_int_params, sym_shared_params, sym_noopt_params, sym_duration});

	% Clean up our symbolic variables
	disp('		Cleaning up interval symbolic variables')
	syms b e i s n t clear
end