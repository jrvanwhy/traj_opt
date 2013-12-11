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

function phase = traj_setup_dynamic_phase(name, dynsys_fcn, state, input, n_intervals, add_params, noopt_params, technique)
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

	% The default value for n_intervals is 50
	if nargin < 5 || isempty(n_intervals)
		n_intervals = 50;
	end

	% The default value for add_params and noopt_params are empty cell arrays
	if nargin < 6 || isempty(add_params)
		add_params = {};
	end
	if nargin < 7 || isempty(noopt_params)
		noopt_params = {};
	end

	% Default to midpoint integration
	if nargin < 8 || isempty(technique)
		technique = 'midpoint';
	end

	% This is the input sanity verification section
	% name should be a character array (string)
	if ~ischar(name)
		error('name must be empty or a string')
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

	% Check that the name is a string
	if ~ischar(technique)
		error('Optimization technique must be a string')
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

	% Set the technique
	phase.technique = technique;

	% Do setup (currently, we only have trapezoidal direct collocation)
	phase = setup_ode_solver(phase, dynsys_fcn);

	% Let the user know when this function terminates
	disp(['Setup for phase ''' name ''' completed.'])
end

% This function calls the dynamic system function and places its outputs into the appropriate
% locations (as symbolic expressions). It needs symbolic samples for the state, input, additional parameters, and non-optimized parameters.
function phase = call_dynsys_fcn(phase, dynsys_fcn, sym_state, sym_input, sym_add_params, sym_noopt_params)
	% Diagnostic message to inform the user of the framework's progress
	disp('		Processing dynamic system function')

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
	disp('			Calling dynamic system function')
	[phase.dynsys.dx,dynsys_varouts{1,1:nargout(dynsys_fcn)-1}] = dynsys_fcn(dynsys_args{:});

	% Create the dx function representation
	disp('			Creating dx function representation')
	phase.dynsys.dx = matlabFunction(simplify(phase.dynsys.dx), 'vars', {sym_state, sym_input, sym_add_params, sym_noopt_params});

	% Initialize dynamic system fields
	phase.dynsys.functions   = {};

	% Initialize temporary cost and constraint variables
	cost = sym(0);
	c    = sym([]);
	ceq  = sym([]);

	% Iterate through the costs and constraints, adding them to phases
	disp('			Adding costs and constraints to dynamic system representation.')
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
				disp(['				Processing constraint ''' output.name ''''])

				% Append it to the appropriate constraints field.
				if output.con_type == '<='
					c = [c; output.fcn(:)];
				else
					ceq = [ceq; output.fcn(:)];
				end

			case 'cost'
				% Progress message so the user doesn't think it's frozen
				disp(['				Processing cost ''' output.name ''''])

				% Add it to the costs field
				cost = cost + output.fcn;

			case 'function'
				% Let the user know things are happening
				disp(['				Processing function ''' output.name ''''])

				% Convert to an anonymous function
				% Append it to the dynamic system functions
				output.fcn = matlabFunction(simplify(output.fcn), 'vars', {sym_state, sym_input, sym_add_params, sym_noopt_params});
				phase.dynsys.functions{end+1} = output;

			otherwise
				% Spit out an error -- this is not an acceptable structure type
				error(['Structure type ''' output.struct_type ''' not a valid return type for the dynamic system function.'])
		end
	end

	% Create cost and constraint anonymous functions
	disp('			Generating anonymous dynamic system cost function')
	phase.dynsys.cost = matlabFunction(cost, 'vars', {sym_state, sym_input, sym_add_params, sym_noopt_params});

	disp('			Generating anonymous dynamic system inequality constraint function')
	phase.dynsys.c    = matlabFunction(c,    'vars', {sym_state, sym_input, sym_add_params, sym_noopt_params});

	disp('			Generating anonymous dynamic system equality constraint function')
	phase.dynsys.ceq  = matlabFunction(ceq,  'vars', {sym_state, sym_input, sym_add_params, sym_noopt_params});
end

% Generate the symbolic parameters for the optimization
function [sym_states,sym_inputs,sym_duration] = gen_params(state_size, state_grid, input_size, input_grid)
	sym_states   = sym('x', [state_size state_grid]);
	sym_inputs   = sym('u', [input_size input_grid]);
	sym_duration = sym('d', 'real');
end

% This function sets up for midpoint-based direct collocation
function phase = setup_dircol_midpoint(phase, dynsys_fcn)
	disp('	Setting up midpoint ODE solver')

	% Set up states, inputs, and duration
	disp('		Creating symbolic parameters')
	[sym_states,sym_inputs,sym_duration] =  ...
		gen_params(numel(phase.names.state), 2, ...
		           numel(phase.names.input), 1);

	% Set up dt
	sym_dt = sym_duration/phase.n_intervals;

	% Generate dynamic system function expressions.
	phase = call_dynsys_fcn(phase, dynsys_fcn, sym_states(:,1), sym_inputs(:,1), sym([]), sym([]));

	% Generate important intermediate values
	mid_states = (sym_states(:,1:end-1) + sym_states(:,2:end))/2;

	% Generate parameter mappings
	disp('		Creating parameter mappings')
	state_params    = zeros(numel(phase.names.state), phase.n_intervals+1);
	input_params    = zeros(numel(phase.names.input), phase.n_intervals);
	state_params(:) = 1:numel(state_params);
	n_params        = numel(state_params);
	input_params(:) = n_params+1:n_params+numel(input_params);
	n_params        = n_params + numel(input_params);
	duration_param  = n_params + 1;
	n_params        = n_params + 1;
	phase.n_params  = n_params;

	intvl_params = [ state_params(:,1:end-1)
	                 state_params(:,2:end)
	                 input_params
	                 duration_param * ones(1, phase.n_intervals) ];
	sym_intvl_params = [ sym_states(:,1)
	                     sym_states(:,2)
	                     sym_inputs
	                     sym_duration ];

	% Generate cost
	phase_cost = sym_dt * phase.dynsys.cost(mid_states, sym_inputs, [], []);
	phase.cost = opt_create_vecfcn(phase_cost, intvl_params, sym_intvl_params);

	disp('		Creating and processing constraints')
	c   = phase.dynsys.c(  mid_states, sym_inputs, [], []);
	ceq = phase.dynsys.ceq(mid_states, sym_inputs, [], []);

	% Generate nonnegative duration constraint
	c = [c; -sym_duration];

	% Create inequality constraint vector function
	phase.c = opt_create_vecfcn(c, intvl_params, sym_intvl_params);

	% Generate collocation constraint
	dstates = phase.dynsys.dx(mid_states, sym_inputs, [], []);
	ceq = [ceq; sym_dt * dstates + sym_states(:,1) - sym_states(:,2)];

	% Create equality constraint vector function
	phase.ceq = opt_create_vecfcn(ceq, intvl_params, sym_intvl_params);

	% Define a few useful functions.
	disp('		Processing functions')
	functions = phase.dynsys.functions;
	functions{end+1} = traj_create_function('dstates', dstates);
	functions{end+1} = traj_create_function('inputs', sym_inputs);

	% Initialize the phases's functions list
	phase.functions = {};

	% Iterate through the functions, creating vecfcn forms for each
	for iterfcn = 1:numel(functions)
		% Convenience variable
		cur_fcn = functions{iterfcn};

		% Evaluate the function expression, if necessary
		if isa(cur_fcn.fcn, 'function_handle')
			cur_fcn.fcn = cur_fcn.fcn(mid_states, sym_inputs, [], []);
		end

		% Create a vecfcn form; make it simple because this isn't used
		% for optimization
		cur_fcn.fcn = opt_create_vecfcn(cur_fcn.fcn, intvl_params, sym_intvl_params, false);

		% Append this function to the phase functions list
		phase.functions{end+1} = cur_fcn;
	end

	% States and duration are special; deal with them separately
	phase.functions{end+1} = traj_create_function('states', sym_states(:,1));
	phase.functions{end}.fcn = opt_create_vecfcn(phase.functions{end}.fcn, ...
		state_params, sym_states(:,1));

	phase.functions{end+1} = traj_create_function('duration', sym_duration);
	phase.functions{end}.fcn = opt_create_vecfcn(phase.functions{end}.fcn, ...
		duration_param, sym_duration);

	% Clean up symbolic variables
	disp('	Cleaning up symbolic variables')
	syms sym_states sym_inputs sym_duration clear
end

% This function determines the appropriate solver to use and calls it
function phase = setup_ode_solver(phase, dynsys_fcn)
	switch phase.technique
		case 'midpoint'
			phase = setup_dircol_midpoint(phase, dynsys_fcn);

		otherwise
			error(['Optimization technique ''' phase.technique ''' not recognized'])
	end
end
