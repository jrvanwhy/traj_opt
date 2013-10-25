% This function prepares the given scenario for optimization. The structure of the scenario struct is given in
% traj_optimize.m

function scenario = traj_setup(scenario)
	disp('Beginning optimization setup')
	scenario = traj_clean_scenario(scenario);
	scenario = traj_gen_system(scenario);
	scenario = traj_gen_sys_bnds(scenario);
	scenario = traj_gen_add_params(scenario);
	% Set the method in the structure (in the future, multiple methods will be available, so there will be an intelligent
	% decision made to choose between the methods).
	scenario.method = 'direct collocation';
	scenario = traj_dircol(scenario);
	scenario = traj_simplify_constraints(scenario);
end

% This cleans out the scenario, so it can be re-setup
function scenario = traj_clean_scenario(scenario)
	scenario.params.params = [];
	scenario.linear_constraints.A   = [];
	scenario.linear_constraints.b   = [];
	scenario.linear_constraints.Aeq = [];
	scenario.linear_constraints.beq = [];
end

% This sets up the cost and collocation constraint functions for direct collocation.
function scenario = traj_dircol(scenario)
	disp('Setting up for Direct Collocation')

	% Generate state, inputs, and duration, and copy them over to the parameters
	scenario.states.states(:,1:2:2*scenario.shoots+1) = sym('x', [size(scenario.states.def_state, 1),   scenario.shoots+1]);
	scenario.inputs.inputs                            = sym('u', [size(scenario.inputs.def_input, 1), 2*scenario.shoots+1]);
	scenario.duration                                 = sym('t', 'real');
	assume(scenario.duration >= 0);
	error(' Old code!')
	state_params = scenario.states.states(:,1:2:2*scenario.shoots+1);
	scenario.params.params = [scenario.params.params; state_params(:); scenario.inputs.inputs(:); scenario.duration];
	scenario.params = declare_field(scenario.params, 'init_params');
	% Only initialize it if it's not already the right size
	if size(scenario.params.init_params, 1) ~= size(scenario.params.params, 1)
		scenario.params.init_params = [repmat(scenario.states.def_state,   scenario.shoots+1, 1)
		                               repmat(scenario.inputs.def_input, 2*scenario.shoots+1, 1)
		                               1];
	end

	% Find the derivative of the state at the grid points
	scenario.states.dstates(:,1:2:2*scenario.shoots+1) = eval_sym_fcn(scenario.system.dx, [scenario.states.states(:,1:2:2*scenario.shoots+1)
	                                                                                       scenario.inputs.inputs(:,1:2:2*scenario.shoots+1)
	                                                                                       repmat(scenario.add_params.params, 1, scenario.shoots+1)]);
	% Evaluate the states and derivative of the states at the midpoints of the grid
	scenario.states.approx_dstates = scenario.states.dstates;
	[scenario.states.states(:,2:2:2*scenario.shoots),scenario.states.approx_dstates(:,2:2:2*scenario.shoots)] = ...
		traj_dircol_mid_approx(scenario.states.states( :,1:2:2*scenario.shoots-1), ...
		                       scenario.states.dstates(:,1:2:2*scenario.shoots-1), ...
		                       scenario.states.states( :,3:2:2*scenario.shoots+1), ...
		                       scenario.states.dstates(:,3:2:2*scenario.shoots+1), scenario.duration / scenario.shoots);

	% Evaluate actual dstates at the midpoints of the shoots
	scenario.states.dstates(:,2:2:2*scenario.shoots) = eval_sym_fcn(scenario.system.dx, [scenario.states.states(:,2:2:2*scenario.shoots)
	                                                                                     scenario.inputs.inputs(:,2:2:2*scenario.shoots)
	                                                                                     repmat(scenario.add_params.params, 1, scenario.shoots)]);

	% Evaluate cost
	scenario.objective.dcost = eval_sym_fcn(scenario.system.dcost, [scenario.states.states
	                                                                scenario.inputs.inputs
	                                                                repmat(scenario.add_params.params, 1, 2*scenario.shoots+1)]);
	scenario.objective.value = 2 * scenario.duration * (sum(scenario.objective.dcost(1:2:2*scenario.shoots+1)) + ...
	                                                2 * sum(scenario.objective.dcost(2:2:2*scenario.shoots))   - ...
	                                               .5 * sum(scenario.objective.dcost([1 2*scenario.shoots+1]))) / scenario.shoots;

	% Evaluate linear inequality and equality constraints
	scenario.constraints.c   = eval_sym_fcn(scenario.system.c,   [scenario.states.states
	                                                              scenario.inputs.inputs
	                                                              repmat(scenario.add_params.params, 1, 2*scenario.shoots+1)]);
	scenario.constraints.ceq = eval_sym_fcn(scenario.system.ceq, [scenario.states.states
	                                                              scenario.inputs.inputs
	                                                              repmat(scenario.add_params.params, 1, 2*scenario.shoots+1)]);
	% We need to rearrange these into column vectors
	scenario.constraints.c   = scenario.constraints.c(:);
	scenario.constraints.ceq = scenario.constraints.ceq(:);

	% Add on collocation constraints
	scenario.constraints.colloc = scenario.states.dstates(:,2:2:2*scenario.shoots) - scenario.states.approx_dstates(:,2:2:2*scenario.shoots);
	scenario.constraints = declare_field(scenario.constraints, 'ceq');
	scenario.constraints.ceq = [scenario.constraints.ceq; scenario.constraints.colloc(:)];

	% Constraint duration to be nonnegative
	scenario.constraints.c = [scenario.constraints.c; 1e-6-scenario.duration];

	% Evaluate the static system function, if defined
	scenario = declare_field(scenario, 'static_sys_fcn');
	if ~isempty(scenario.static_sys_fcn)
		[static_cost,static_c,static_ceq] = scenario.static_sys_fcn(scenario);
		scenario.objective.value = scenario.objective.value + static_cost;
		scenario.constraints.c   = [scenario.constraints.c;   static_c];
		scenario.constraints.ceq = [scenario.constraints.ceq; static_ceq];
	end

	% Set up the jacobians (we'll ignore the hessians for now, since they're not used and are complex to compute).
	scenario.objective.jacobian  = jacobian(scenario.objective.value, scenario.params.params);
	scenario.constraints.c_jac   = jacobian(scenario.constraints.c,   scenario.params.params);
	scenario.constraints.ceq_jac = jacobian(scenario.constraints.ceq, scenario.params.params);
end

% This function evaluates the midpoint approximation for direct collocation
function [x_mid,dx_mid] = traj_dircol_mid_approx(x0, dx0, x1, dx1, dt)
	% Reshape the inputs, placing them into a matrix, then use linear algebra to get the result. Lastly, reshape the outputs into the original shape
	x_mid  = reshape([x0(:), dx0(:), x1(:), dx1(:)] * [     .5; dt/8;     .5; dt/-8], size(x0));
	dx_mid = reshape([x0(:), dx0(:), x1(:), dx1(:)] * [-1.5/dt; -.25; 1.5/dt;  -.25], size(x0));
end

% This function generates the symbolic additional parameters
function scenario = traj_gen_add_params(scenario)
	scenario = declare_field(scenario, 'add_params');
	scenario.add_params = declare_field(scenario.add_params, 'value');
	scenario.add_params.params = sym('p', size(scenario.add_params.value));
end

% This function evaluates the given symbolic function with the given arguments.
% The arguments should be specified as a matrix -- each column calls the function once, and each row represents one value of one
% of the function's arguments.
% If fcn or inputs is empty, this returns an empty array.
function out = eval_sym_fcn(fcn, inputs)
	if isempty(fcn) || isempty(inputs)
		out = [];
		return;
	end

	% Found at http://www.mathworks.com/matlabcentral/answers/85539-creating-and-differentiating-symbolic-functions-with-vector-inputs-of-different-dimensions
	% There really should be a better way...
	input_cell = num2cell(inputs);
	for i = 1:size(inputs,2)
		out(:,i) = fcn(input_cell{:,i});
	end
end

% This function checks if the dimensions of the system results make sense.
% If it catches an issue, it throws an error
function traj_check_system(scenario, dx, dcost, c, ceq)
	% Check dx
	if size(dx) ~= size(scenario.states.def_state)
		error('System dx has an incorrect size.')
	end

	% Check dcost
	if size(dcost) ~= [1 1]
		error('System dcost is not a scalar.')
	end

	% Check c
	if ~ismember(size(c, 2), [0 1])
		error('System c is not a column vector.')
	end

	% Check ceq
	if ~ismember(size(ceq, 2), [0 1])
		error('System ceq is not a column vector.')
	end
end

% This function generates the bounds for the inputs and states, if not user-specified
function scenario = traj_gen_sys_bnds(scenario)
	disp('Generating default system bounds')

	scenario.inputs = declare_field(scenario.inputs, 'lb', -Inf * ones(size(scenario.inputs.def_input)));
	scenario.inputs = declare_field(scenario.inputs, 'ub',  Inf * ones(size(scenario.inputs.def_input)));
	scenario.states = declare_field(scenario.states, 'lb', -Inf * ones(size(scenario.states.def_state)));
	scenario.states = declare_field(scenario.states, 'ub',  Inf * ones(size(scenario.states.def_state)));
end

% This function generates the symbolic system function from the actual system function
function scenario = traj_gen_system(scenario)
	disp('Generating symbolic system function')

	% Generate our inputs
	state      = sym('x', size(scenario.states.def_state));
	input      = sym('u', size(scenario.inputs.def_input));
	add_params = sym('p', size(scenario.add_params.value));

	% Call the system function, to get symbolic representations for each of its outputs
	[dx,dcost,c,ceq] = scenario.system_fcn(state, input, add_params);

	% Check if the system's output makes sense (dimensionwise)
	traj_check_system(scenario, dx, dcost, c, ceq);

	% Convert the symbolic representations into symbolic functions
	scenario.system.dx    = symfun(dx,    [state; input; add_params]);
	scenario.system.dcost = symfun(dcost, [state; input; add_params]);
	% c and ceq may be empty -- if so, don't generate them
	scenario.system.c   = [];
	scenario.system.ceq = [];
	if ~isempty(c)
		scenario.system.c     = symfun(c,     [state; input; add_params]);
	end
	if ~isempty(ceq)
		scenario.system.ceq   = symfun(ceq,   [state; input; add_params]);
	end

	% Clean up our symbolic variables before they leave the workspace
	syms state input add_params clear
end

% This function simplifies the constraint functions. It locates linear constraints and moves them to the
% linear constraints struct in scenario. It then sorts through this struct, identifying specific forms of
% linear constraints, doing specific optimizations for each form.
function scenario = traj_simplify_constraints(scenario)
	disp('Simplifying constraints')
	scenario                    = declare_field(scenario, 'linear_constraints');
	scenario.linear_constraints = declare_field(scenario.linear_constraints, 'A');
	scenario.linear_constraints = declare_field(scenario.linear_constraints, 'b');
	scenario.linear_constraints = declare_field(scenario.linear_constraints, 'Aeq');
	scenario.linear_constraints = declare_field(scenario.linear_constraints, 'beq');
	done = false;
	while (~done)
		c_lin_cons   = traj_id_lin_cons([scenario.constraints.c_jac],   scenario.add_params.params);
		ceq_lin_cons = traj_id_lin_cons([scenario.constraints.ceq_jac], scenario.add_params.params);

		disp(['	Found ', num2str(sum(c_lin_cons) + sum(ceq_lin_cons)), ' linear constraints out of ', ...
			num2str(size(c_lin_cons, 1) + size(ceq_lin_cons, 1)), ' constraints total.'])

		[scenario.linear_constraints.A,scenario.linear_constraints.b] = ...
			traj_conv_lin_cons(scenario.linear_constraints.A, ...
			                   scenario.linear_constraints.b, ...
			                   scenario.constraints.c,        ...
			                   scenario.constraints.c_jac,    ...
			                   c_lin_cons,                    ...
			                   scenario);

		[scenario.linear_constraints.Aeq,scenario.linear_constraints.beq] = ...
			traj_conv_lin_cons(scenario.linear_constraints.Aeq, ...
			                   scenario.linear_constraints.beq, ...
			                   scenario.constraints.ceq,        ...
			                   scenario.constraints.ceq_jac,    ...
			                   ceq_lin_cons,                    ...
			                   scenario);

		[scenario.constraints.c,scenario.constraints.c_jac] = ...
			traj_clean_lin_cons(scenario.constraints.c, scenario.constraints.c_jac, c_lin_cons);
		[scenario.constraints.ceq,scenario.constraints.ceq_jac] = ...
			traj_clean_lin_cons(scenario.constraints.ceq, scenario.constraints.ceq_jac, ceq_lin_cons);

		% Identify and remove fixed params.
		[scenario,num_fixed_params] = traj_remove_fixed_params(scenario);

		% Remove empty constraints that may have been generated as a result of fixed
		% parameter substitution.
		scenario = traj_remove_empty_cons(scenario);

		% If we found any fixed params, we should re-run, since there may be new linear constraints and fixed params
		if num_fixed_params
			disp(['	Since there were fixed parameters, re-running constraint simplification to catch ', ...
				'new simplifiable constraints.'])
			done = false;
		else
			done = true;
		end
	end
end

% This function identifies linear constraints, given the jacobian of a set of constraints
% It returns a vector of 1's and 0's -- 1 means linear, 0 means nonlinear. This works for both
% equality and inequality constraints. It works by looking at the rows of the jacobian -- if a row
% is constant with respect to the optimization parameters, then it is linear. Note that the params
% input refers to the additional parameter list, scenario.add_params.params
function lin_constraints = traj_id_lin_cons(jacobian, add_params)
	% symvar is not vectorizable, therefore we need to iterate through the constraints
	parfor i = 1:size(jacobian, 1)
		% Check if the given constraint is linear by checking if its jacobian is only a function of the
		% additional parameters.
		vars = symvar(jacobian(i,:));
		lin_constraints(i,1) = isempty(setdiff(vars, add_params));
	end
end

% This function removes the given linear constraints from the nonlinear constraint function.
% con is the nonlinear constraint function, con_jac is its jacobian, and lin_cons is the list of linear constraints
% returned by traj_id_lin_cons
function [con,con_jac] = traj_clean_lin_cons(con, con_jac, lin_cons)
	% Locate nonlinear constraints and keep those
	nonlin_cons = find(~lin_cons);

	% Copy over the nonlinear constraints
	con     = con(nonlin_cons,:);
	con_jac = con_jac(nonlin_cons,:);
end

% This function turns linear "nonlinear" constraints into equivalent linear constraints
% A and b are the existing linear constraints (they will be appended to).
% con is the constraint function; con_jac is its jacobian
% lin_cons is the vector of linear and nonlinear constraints, as returned by traj_id_lin_cons
% This returns the A and b functions for the linear constraints
function [A,b] = traj_conv_lin_cons(A, b, con, con_jac, lin_cons, scenario)
	% Identify the indices of the linear constraints
	lin_con_inds = find(lin_cons);

	% The A matrix is just the relevant rows of the jacobian of the constraints
	A = [A; con_jac(lin_con_inds,:)];

	% The b matrix is the negative of the relevant rows of the constraint function evaluated at zero
	b = [b; subs(-con(lin_con_inds), scenario.params.params, zeros(size(scenario.params.params)))];
end

% This function removes empty linear constraints from the scenario. Empty linear constraints
% can occur as a result of the fixed parameter substitution process. Empty linear constraints
% have all 0 elements in their row in the A matrix
function scenario = traj_remove_empty_cons(scenario)
	disp('	Removing empty constraints that may have been produced as a result of fixed parameter substitution.')

	% This function identifies nonempty constraints, returning a vector with row indices for each nonempty constraint.
	function indices = id_useful_cons(A)
		% The sum errors out if A is empty, so check for that condition
		if isempty(A)
			indices = [];
		else
			indices = find(sum(A, 2));
		end
	end

	% Clean out inequality and equality constraints separately (but in parallel).
	c_useful   = id_useful_cons(scenario.linear_constraints.A);
	ceq_useful = id_useful_cons(scenario.linear_constraints.Aeq);

	% ... and the actual cleanup.
	scenario.linear_constraints.A   = scenario.linear_constraints.A(  c_useful,   :);
	scenario.linear_constraints.b   = scenario.linear_constraints.b(  c_useful,   :);
	scenario.linear_constraints.Aeq = scenario.linear_constraints.Aeq(ceq_useful, :);
	scenario.linear_constraints.beq = scenario.linear_constraints.beq(ceq_useful, :);
end

% This function removes parameters that were fixed to some value by a constraint.
% It returns fixed_params so that higher level code can re-run the constraint simplification
function [scenario, num_fixed_params] = traj_remove_fixed_params(scenario)
	disp('	Identifying and removing fixed parameters.');

	% Identify the constraints corresponding to the fixed params (and set up variables to be used later).
	A_logicals = ~logical(scenario.linear_constraints.Aeq == 0);
	fixed_cons_A_raw = sum(A_logicals, 2) == 1;
	fixed_cons_A_rows = find(fixed_cons_A_raw);
	fixed_cons_A_log = A_logicals(fixed_cons_A_rows,:);

	% Identify the indices of the fixed params
	[fixed_con_indxs,fixed_param_indxs] = find(fixed_cons_A_log);

	% Let the user know if this worked.
	num_fixed_params = size(fixed_param_indxs, 1);
	disp(['	Identified ', num2str(num_fixed_params), ' fixed parameters.'])

	% Abort if there are 0 fixed params
	if ~num_fixed_params
		scenario = scenario;
		return;
	end

	% Find the fixed parameters themselves
	fixed_params = scenario.params.params(fixed_param_indxs);

	% Now compute the values of the fixed parameters
	fixed_param_values = scenario.linear_constraints.beq(fixed_cons_A_rows(fixed_con_indxs)) ...
		./ scenario.linear_constraints.Aeq( ...
		sub2ind(size(scenario.linear_constraints.Aeq), fixed_cons_A_rows(fixed_con_indxs), fixed_param_indxs));

	% A few more variables we'll need later
	fixed_param_vec = sum(fixed_cons_A_log, 1);
	fixed_param_idxs    = find(fixed_param_vec);
	nonfixed_param_idxs = find(~fixed_param_vec);

	% Shrink the linear constraints, removing fixed parameters and their constraints
	nonfixed_cons_indxs = find(~fixed_cons_A_raw);

	[scenario.linear_constraints.A,scenario.linear_constraints.b] = ...
		traj_fixed_param_lin_con(scenario.linear_constraints.A, ...
		                         scenario.linear_constraints.b, ...
		                         fixed_param_vec,               ...
		                         fixed_param_values);

	[scenario.linear_constraints.Aeq,scenario.linear_constraints.beq] = ...
		traj_fixed_param_lin_con(scenario.linear_constraints.Aeq, ...
		                         scenario.linear_constraints.beq, ...
		                         fixed_param_vec,                 ...
		                         fixed_param_values,              ...
		                         nonfixed_cons_indxs);

	% Replace the parameters everywhere else they're found.
	fprintf('	Done computing parameter replacements, replacing fixed parameters')
	scenario = traj_replace_fixed_params(scenario, fixed_params, fixed_param_values, nonfixed_param_idxs);
	disp('	Done replacing fixed parameters.')

	disp('	Removing fixed parameters from parameters list.')
	sym(scenario.params.params(fixed_param_idxs), 'clear');
	scenario.params.params = scenario.params.params(nonfixed_param_idxs);
	scenario.params.init_params = scenario.params.init_params(nonfixed_param_idxs);
end

% This function replaces the fixed parameters in the various expressions in scenario
function scenario = traj_replace_fixed_params(scenario, fixed_params, fixed_param_values, nonfixed_param_idxs)
	% Execute the replacements in parallel (I wish there was a better way...)
	% We pass in this struct with information on the fixed params; this is just for convenience
	paramvals.fixed_params        = fixed_params;
	paramvals.fixed_param_values  = fixed_param_values;
	paramvals.nonfixed_param_idxs = nonfixed_param_idxs;
	% Pull everything out of scenario; else the parfor loop won't compile
	objective_jacobian    = scenario.objective.jacobian;
	constraints_ceq_jac   = scenario.constraints.ceq_jac;
	constraints_c_jac     = scenario.constraints.c_jac;
	constraints_ceq       = scenario.constraints.ceq;
	states_states         = scenario.states.states;
	states_dstates        = scenario.states.dstates;
	objective_dcost       = scenario.objective.dcost;
	states_approx_dstates = scenario.states.approx_dstates;
	constraints_c         = scenario.constraints.c;
	inputs_inputs         = scenario.inputs.inputs;
	objective_value       = scenario.objective.value;
	duration              = scenario.duration;
	for i = 1:12
		% Pretty dots for the user
		fprintf('.')

		% These are sorted from slowest to fastest to maximize parallel processor utilization
		switch i
			case  1; objective_jacobian    = traj_fixed_params_fix_jacobian(paramvals, objective_jacobian);
			case  2; constraints_ceq_jac   = traj_fixed_params_fix_jacobian(paramvals, constraints_ceq_jac);
			case  3; constraints_c_jac     = traj_fixed_params_fix_jacobian(paramvals, constraints_c_jac);
			case  4; constraints_ceq       = traj_expr_replace_fixed_params(paramvals, constraints_ceq);
			case  5; states_states         = traj_expr_replace_fixed_params(paramvals, states_states);
			case  6; states_dstates        = traj_expr_replace_fixed_params(paramvals, states_dstates);
			case  7; objective_dcost       = traj_expr_replace_fixed_params(paramvals, objective_dcost);
			case  8; states_approx_dstates = traj_expr_replace_fixed_params(paramvals, states_approx_dstates);
			case  9; constraints_c         = traj_expr_replace_fixed_params(paramvals, constraints_c);
			case 10; inputs_inputs         = traj_expr_replace_fixed_params(paramvals, inputs_inputs);
			case 11; objective_value       = traj_expr_replace_fixed_params(paramvals, objective_value);
			case 12; duration              = traj_expr_replace_fixed_params(paramvals, duration);
		end
	end
	% Now we stuff it all back into scenario
	scenario.objective.jacobian    = objective_jacobian;
	scenario.constraints.ceq_jac   = constraints_ceq_jac;
	scenario.constraints.c_jac     = constraints_c_jac;
	scenario.constraints.ceq       = constraints_ceq;
	scenario.states.states         = states_states;
	scenario.states.dstates        = states_dstates;
	scenario.objective.dcost       = objective_dcost;
	scenario.states.approx_dstates = states_approx_dstates;
	scenario.constraints.c         = constraints_c;
	scenario.inputs.inputs         = inputs_inputs;
	scenario.objective.value       = objective_value;
	scenario.duration              = duration;
	fprintf('\n')
end

% This function replaces the fixed parameters in the given expression, then returns it
function expr = traj_expr_replace_fixed_params(paramvals, expr)
	% For whatever reason, substituting them all at once takes forever (and arbitrary memory). Substitute one at a time
	% instead for performance.
	for j = 1:size(paramvals.fixed_params, 1)
		% According to the matlab documentation, if the "old" value (paramvals.fixed_params(j)) is not
		% in expr, it'll switch its parameters around (‽)
		% Check to make sure that this value's actually in expr, or we get weird random stuff.
		if ~ismember(paramvals.fixed_params(j), symvar(expr))
			continue
		end

		expr = subs(expr, paramvals.fixed_params(j), paramvals.fixed_param_values(j));
	end
end

% This function fixes a jacobian
function jac = traj_fixed_params_fix_jacobian(paramvals, jac)
	jac = traj_expr_replace_fixed_params(paramvals, jac(:,paramvals.nonfixed_param_idxs));
end

% This removes fixed parameters from linear constraints
% A and b are the linear constraint matrices
% fixed_param_vec is a vector of 1's and 0's, identifying fixed (1) and non-fixed (0) parameters
% fixed_param_vals is the vector of values of the fixed parameters
% nonfixed_cons_indxs should only be specified for the equality constraints -- it indicates
% which constraints correspond to non-fixed-param constraints
function [A,b] = traj_fixed_param_lin_con(A, b, fixed_param_vec, fixed_param_vals, nonfixed_cons_indxs)
	if nargin >= 5
		% Remove constraints corresponding to the fixed parameters
		A = A(nonfixed_cons_indxs,:);
		b = b(nonfixed_cons_indxs);
	end

	% Update b to reflect the removal of the fixed parameters
	b = b - A(:,find(fixed_param_vec)) * fixed_param_vals;

	% Remove the fixed parameters from A
	A = A(:, find(~fixed_param_vec));
end
