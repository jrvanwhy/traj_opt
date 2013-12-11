% This function does most of the setup and processing for a given scenario.
% It returns a scenario ready for optimization (or generation
% of a fast optimization function, if desired).

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

	% Call the scenario function and process its outputs
	scenario = process_scenario_fcn(scenario, scenario_fcn);

	% Combine all the optimization functions into one big list
	scenario = combine_optfcns(scenario);
end

% This function adds a phase to the given scenario
function scenario = add_phase(scenario, phase)
	% Append to the phases list (creating it if necessary)
	if isfield(scenario, 'phases')
		scenario.phases{end+1} = phase;
	else
		scenario.phases{1} = phase;
	end
end

% This function combines all optimization functions into one big list for optimization
function scenario = combine_optfcns(scenario)
	% Initialize the output arrays to the scenario function outputs
	scenario.costs = {scenario.scenfun.cost};
	scenario.c     = {scenario.scenfun.c};
	scenario.ceq   = {scenario.scenfun.ceq};

	% Iterate through the phases, copying the functions from each phase
	for iterphase = 1:numel(scenario.phases)
		phase = scenario.phases{iterphase};
		disp(['	Combining functions from phase ''' phase.names.phase ''''])

		% Copy over costs, c, and ceq
		scenario.costs{end+1} = phase.cost;
		scenario.c{end+1}     = phase.c;
		scenario.ceq{end+1}   = phase.ceq;
	end
end

% This function finds variables within the optimization parameters
function [locs,vars] = find_vars(vars, opt_params)
	[vars,~,locs] = intersect(vars, opt_params, 'R2012a');
end

% This function generates the initial guess
function x0 = gen_x0(init_params, init_vals, opt_params)
	x0 = zeros(numel(opt_params), 1);

	% Iterate through the initial values, setting the relevant
	% parameters
	for iterinit = 1:numel(init_params)
		param = init_params(iterinit);
		value = init_vals(iterinit);
		x0(logical(param == opt_params)) = value;
	end
end

% Sets up parameter maps for the phases
function scenario = process_phases(scenario)
	cur_params = 0;
	% Iterate through parameter counts, putting ranges into scenario.param_map
	for iterphase = 1:numel(scenario.phases)
		phase = scenario.phases{iterphase};
		disp(['	Processing parameters for phase ''' phase.names.phase ''''])
		scenario.param_map(1, iterphase) = cur_params + 1;
		scenario.param_map(2, iterphase) = cur_params + phase.n_params;
		cur_params = cur_params + phase.n_params;
	end

	% Save the number of parameters
	scenario.n_params = cur_params;
end

% Handles everything necessary for the scenario function
function scenario = process_scenario_fcn(scenario, scenario_fcn)
	% Check if a scenario function was given. If not, return scenario.
	if isempty(scenario_fcn)
		scenario = scenario;
		return
	end

	disp('	Setting up to call scenario function')

	% We need to generate the full optimization parameters list
	% so we can symbolically evaluate all functions
	disp('	Generating optimization parameters')
	opt_params = sym('x', [scenario.n_params 1]);

	% Add the symbolic outputs for each function to the scenario
	scenario = traj_eval_funcs(scenario, opt_params, []);

	% Call the scenario function, capturing all output arguments
	disp('	Calling scenario function')
	[scenario_varouts{1:nargout(scenario_fcn)}] = scenario_fcn(scenario);

	% Initialize outputs
	cost        = sym(0);
	c           = sym([]);
	ceq         = sym([]);
	init_params = sym([]);
	init_vals   = [];

	% Set initial durations to 1, to avoid infeasibility of the linearized parameters
	for iterphase = 1:numel(scenario.phases)
		phase = scenario.phases{iterphase};
		init_params(end+1) = phase.duration;
		init_vals(end+1)   = 1;
	end

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

				% Append the constraint based on its type
				if strcmp(output.con_type, '<=')
					c = [c; output.fcn(:)];
				else
					ceq = [ceq; output.fcn(:)];
				end

			case 'cost'
				% Progress message so the user doesn't think it's frozen
				disp(['		Processing cost ''' output.name ''''])

				% Add this to the symbolic cumulative cost
				cost = cost + sum(output.fcn(:));

			case 'initval'
				% Simply append its parameter and value to the initial values lists
				init_params = [init_params; output.params];
				init_vals   = [init_vals;   output.values];

			otherwise
				% Spit out an error -- this is not an acceptable structure type
				error(['Structure type ''' output.struct_type ''' not a valid return type for the dynamic system function.'])
		end
	end

	% Generate the cost and constraint vecfcns
	disp('	Creating scenario function cost vecfcn')
	[cost_params,sym_cost_params] = find_vars(symvar(cost), opt_params);
	scenario.scenfun.cost = opt_create_vecfcn(cost, cost_params, sym_cost_params);

	disp('	Creating scenario function inequality constraint vecfcn')
	[c_params,sym_c_params] = find_vars(symvar(c), opt_params);
	scenario.scenfun.c = opt_create_vecfcn(c, c_params, sym_c_params);

	disp('	Creating scenario function equality constraint vecfcn')
	[ceq_params,sym_ceq_params] = find_vars(symvar(ceq), opt_params);
	scenario.scenfun.ceq = opt_create_vecfcn(ceq, ceq_params, sym_ceq_params);

	% Create the initial guess using the user's supplied initial values
	disp('	Creating initial guess')
	scenario.x0 = gen_x0(init_params, init_vals, opt_params);

	disp(['	Cleaning up symbolic variables'])
	syms opt_params clear
end
