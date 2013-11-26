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
	phase.dynsys.constraints = {};
	phase.dynsys.costs       = {};
	phase.dynsys.functions   = {};

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

				% Convert the constraint function to an anonymous function, then append it to the
				% constraints field.
				output.fcn = matlabFunction(simplify(output.fcn), 'vars', {sym_state, sym_input, sym_add_params, sym_noopt_params});
				phase.dynsys.constraints = [phase.dynsys.constraints {output}];

			case 'cost'
				% Progress message so the user doesn't think it's frozen
				disp(['				Processing cost ''' output.name ''''])

				% Convert the cost function to an anonymous function with appropriate inputs, then
				% append it to the costs field
				output.fcn = matlabFunction(simplify(output.fcn), 'vars', {sym_state, sym_input, sym_add_params, sym_noopt_params});
				phase.dynsys.costs{end+1} = output;

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
end

% Do final hessian summation and conversion into an anonymous function
function hess = finalize_hessian(opt_params, hcost, hcons)
	disp('		Finalizing Hessian')
	hess.i  = [hcost.i; hcons.i];
	hess.j  = [hcost.j; hcons.j];
	hess.mn =  hcost.mn;

	hess.s  = matlabFunction([hcost.s; hcons.s], 'vars', {opt_params, hcons.lambdas});

	disp('		Cleaning up symbolic lambdas')
	syms l clear
end

% Generate the symbolic parameters for the optimization
function [sym_states,sym_inputs,sym_duration,opt_params,n_params] = gen_params(state_size, state_grid, input_size, input_grid)
	sym_states   = sym('x', [state_size state_grid]);
	sym_inputs   = sym('u', [input_size input_grid]);
	sym_duration = sym('d', 'real');
	opt_params   = [sym_states(:); sym_inputs(:); sym_duration];
	n_params     = numel(opt_params);
end

% This function processes constraints. It also adds in the nonnegative duration constraint
function [phase,hcons] = process_constraints(phase, opt_params, sym_duration, colloc_con, con_eval_fcn)
	% Just copy over the collocation constraint, since it's symbolic already
	disp('		Handling collocation constraint')
	ceq_expr = colloc_con.fcn;

	% Create nonnegative duration constraint
	disp('		Creating nonnegative duration constraint')
	c_expr = -sym_duration;

	% Handle the rest of the constraints
	disp('		Processing user-specified phase constraint functions')
	for itercon = 1:numel(phase.dynsys.constraints)
		con = phase.dynsys.constraints{itercon};
		disp(['			Processing constraint ''' con.name ''''])

		% Detect constraint type
		con_expr = con_eval_fcn(con.fcn);
		if con.con_type == '<='
			c_expr   = [c_expr; con_expr(:)];
		else
			ceq_expr = [ceq_expr; con_expr(:)];
		end
	end
	disp('		Converting inequality constraints to a function')
	phase.c   = matlabFunction(c_expr,   'vars', {opt_params});
	disp('		Converting equality constraints to a function')
	phase.ceq = matlabFunction(ceq_expr, 'vars', {opt_params});
	disp('		Generating inequality constraint jacobian expression')
	jc_expr   = jacobian(c_expr,   opt_params);
	disp('		Generating equality constraint jacobian expression')
	jceq_expr = jacobian(ceq_expr, opt_params);
	disp('		Converting inequality constraint jacobian into sparse form')
	[phase.jc.m,phase.jc.n]                 = size(jc_expr);
	[phase.jc.i,phase.jc.j,jc_s_expr]       = find(jc_expr);
	disp('		Converting equality constraint jacobian into sparse form')
	[phase.jceq.m,phase.jceq.n]             = size(jceq_expr);
	[phase.jceq.i,phase.jceq.j,jceq_s_expr] = find(jceq_expr);
	disp('		Creating inequality constraint jacobian value function')
	phase.jc.s                        = matlabFunction(jc_s_expr,   'vars', {opt_params});
	disp('		Creating equality constraint jacobian value function')
	phase.jceq.s                      = matlabFunction(jceq_s_expr, 'vars', {opt_params});

	disp('		Generating contraint hessian')
	disp('			Creating symbolic lambdas')
	n_c     = numel(c_expr);
	n_ceq   = numel(ceq_expr);
	n_cons  = n_c + n_ceq;
	lambdas = sym('l', [n_cons 1]);

	disp('			Creating sparse symbolic hessian')
	hcons = traj_gen_sparse_hess(opt_params, lambdas, [jc_expr; jceq_expr]);

	% Return the lambdas too
	hcons.lambdas = lambdas;
end

% This function processes the solver-defined and user-defined functions
function functions = process_functions(functions, opt_params, fcn_eval_fcn)
	% Handle the dynamic system functions (plus our own ones above)
	for iterfcn = 1:numel(functions)
		fcn = functions{iterfcn};
		disp(['			Processing function ''' fcn.name ''''])

		% If it's an anonymous function, then call it to obtain the output expression
		if isa(fcn.fcn, 'function_handle')
			fcn.fcn = fcn_eval_fcn(fcn.fcn);
		end

		% Make it an anonymous function
		fcn.fcn = matlabFunction(fcn.fcn, 'vars', {opt_params});

		% Copy function over
		functions{iterfcn}.fcn = fcn.fcn;
	end
end

% This function sets up for midpoint-based direct collocation
function phase = setup_dircol_midpoint(phase, dynsys_fcn)
	disp('	Setting up midpoint ODE solver')

	% Set up states, inputs, and duration
	disp('		Creating symbolic parameters')
	[sym_states,sym_inputs,sym_duration,opt_params,phase.n_params] =  ...
		gen_params(numel(phase.names.state), phase.n_intervals+1, ...
		           numel(phase.names.input), phase.n_intervals);

	% Set up dt
	sym_dt = sym_duration/phase.n_intervals;

	% Generate dynamic system function expressions.
	phase = call_dynsys_fcn(phase, dynsys_fcn, sym_states(:,1), sym_inputs(:,1), sym([]), sym([]));

	% Generate important intermediate values
	mid_states = (sym_states(:,1:end-1) + sym_states(:,2:end))/2;

	% Generate costs
	function dcost_expr = mid_cost_fcn(cost_fcn)
		% Midpoint integration!
		dcost_expr = sym_dt * cost_fcn(mid_states, sym_inputs, [], []);
	end
	[phase.cost,hcost] = sum_costs(phase.dynsys.costs, opt_params, @mid_cost_fcn);

	% Generate collocation constraint
	disp('		Creating and processing constraints')
	dstates = phase.dynsys.dx(mid_states, sym_inputs, [], []);
	colloc_con = traj_create_constraint('collocation', ...
		sym_dt * dstates, '=', sym_states(:,2:end) - sym_states(:,1:end-1));

	function newcon_expr = mid_con_eval_fcn(con_fcn)
		newcon_expr = con_fcn(mid_states, sym_inputs, [], []);
	end
	[phase,hcons] = process_constraints(phase, opt_params, sym_duration, colloc_con, @mid_con_eval_fcn);

	disp('		Combining cost and constraint hessians')
	phase.hess = finalize_hessian(opt_params, hcost, hcons);

	% Define a few useful functions.
	% For simplicity, we inject them directly into the outputs from the dynamic system function
	disp('		Processing functions')
	phase.functions = phase.dynsys.functions;
	phase.functions{end+1} = traj_create_function('states', sym_states);
	phase.functions{end+1} = traj_create_function('dstates', dstates);
	phase.functions{end+1} = traj_create_function('inputs', sym_inputs);
	phase.functions{end+1} = traj_create_function('duration', sym_duration);

	% Function evaluation function
	function fcn_vals = mid_fcn_eval_fcn(fcn)
		fcn_vals = fcn(mid_states, sym_inputs, [], []);
	end
	phase.functions = process_functions(phase.functions, opt_params, @mid_fcn_eval_fcn);

	% Clean up symbolic variables
	disp('	Cleaning up symbolic variables')
	syms sym_states sym_inputs sym_duration clear
end

% This function sets up for trapezoidal direct collocation
function phase = setup_dircol_trapz(phase, dynsys_fcn)
	disp('	Setting up trapezoidal ODE solver')

	% Set up states, inputs, and duration
	disp('		Creating symbolic parameters')
	[sym_states,sym_inputs,sym_duration,opt_params,phase.n_params] =  ...
		gen_params(numel(phase.names.state), phase.n_intervals+1, ...
		           numel(phase.names.input), phase.n_intervals+1);

	% Set up dt
	sym_dt = sym_duration/phase.n_intervals;

	% Generate dynamic system function expressions.
	phase = call_dynsys_fcn(phase, dynsys_fcn, sym_states(:,1), sym_inputs(:,1), sym([]), sym([]));

	% Generate important intermediate values
	dstates    = phase.dynsys.dx(sym_states, sym_inputs, [], []);
	mid_states = (sym_states(:,1:end-1) + sym_states(:,2:end))/2;
	mid_inputs = (sym_inputs(:,1:end-1) + sym_inputs(:,2:end))/2;

	% Generate costs
	function dcost_expr = trapz_cost_fcn(cost_fcn)
		% Use Simpson's rule integration to eliminate false minima resulting from the discretization.
		dcost_expr = sym_dt * trapz(cost_fcn(sym_states, sym_inputs, [], []))/3 + ...
		             sym_dt * sum(cost_fcn(mid_states, mid_inputs, [], [])) * 2/3;
	end
	[phase.cost,hcost] = sum_costs(phase.dynsys.costs, opt_params, @trapz_cost_fcn);

	% Generate collocation constraint
	disp('		Creating and processing constraints')
	colloc_con = traj_create_constraint('collocation', ...
		sym_dt * (dstates(:,1:end-1) + dstates(:,2:end)), '=', ...
		2 * (sym_states(:,2:end) - sym_states(:,1:end-1)));

	function newcon_expr = trapz_con_eval_fcn(con_fcn)
		newcon_expr = con_fcn(sym_states, sym_inputs, [], []);
	end
	[phase,hcons] = process_constraints(phase, opt_params, sym_duration, colloc_con, @trapz_con_eval_fcn);

	disp('		Combining cost and constraint hessians')
	phase.hess = finalize_hessian(opt_params, hcost, hcons);

	% Define a few useful functions.
	% For simplicity, we inject them directly into the outputs from the dynamic system function
	disp('		Processing functions')
	phase.functions = phase.dynsys.functions;
	phase.functions{end+1} = traj_create_function('states', sym_states);
	phase.functions{end+1} = traj_create_function('inputs', sym_inputs);
	phase.functions{end+1} = traj_create_function('dstates', dstates);
	phase.functions{end+1} = traj_create_function('duration', sym_duration);

	% Function evaluation function
	function fcn_vals = trapz_fcn_eval_fcn(fcn)
		fcn_vals = fcn(sym_states, sym_inputs, [], []);
	end
	phase.functions = process_functions(phase.functions, opt_params, @trapz_fcn_eval_fcn);

	% Clean up symbolic variables
	disp('	Cleaning up symbolic variables')
	syms sym_states sym_inputs sym_duration clear
end

% This function determines the appropriate solver to use and calls it
function phase = setup_ode_solver(phase, dynsys_fcn)
	switch phase.technique
		case 'trapezoidal'
			phase = setup_dircol_trapz(phase, dynsys_fcn);

		case 'midpoint'
			phase = setup_dircol_midpoint(phase, dynsys_fcn);

		otherwise
			error(['Optimization technique ''' phase.technique ''' not recognized'])
	end
end

% This function iterates through and sums up the costs.
% You need to pass in a cost function, since this varies by solver type
function [cost_fcn,hcost] = sum_costs(costs, opt_params, solver_cost_fcn)
	disp('		Creating phase cost expression')
	cost_expr = sym(0);
	for itercost = 1:numel(costs)
		cost = costs{itercost};
		disp(['			Processing cost ''' cost.name ''''])

		cost_expr = cost_expr + sum(solver_cost_fcn(cost.fcn));
	end
	disp('		Converting phase cost expression to a function')
	cost_fcn = matlabFunction(cost_expr, 'vars', {opt_params});

	% Go ahead and generate the hessian of the cost
	disp('		Generating sparse hessian of phase cost')
	hcost = traj_gen_sparse_hess(opt_params, 1, jacobian(cost_expr, opt_params));
end
