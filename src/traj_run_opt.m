% This function runs the optimization for the given scenario. Note that the optimization must have already been set up.

function scenario = traj_run_opt(scenario)
	scenario = traj_gen_fmincon_fcns(scenario);
	scenario = traj_opt_fmincon(scenario);
end

% This generates the functions called by fmincon
function scenario = traj_gen_fmincon_fcns(scenario)
	% Check if we need to generate the functions first
	if (isfield(scenario,           'fmincon') && ...
	    isfield(scenario.fmincon, 'objective') && ...
	    isfield(scenario.fmincon,   'nonlcon'))
		% Make MATLAB happy
		scenario = scenario;
		return;
	end

	% Set up the functions fmincon will be calling
	disp('Generating MATLAB objective function')
	scenario.fmincon.objective = matlabFunction(scenario.objective.value, scenario.objective.jacobian.', 'vars', {scenario.params.params, scenario.add_params.params});

	disp('Generating MATLAB constraint function')
	scenario.fmincon.nonlcon   = matlabFunction(scenario.constraints.c,       scenario.constraints.ceq, ...
	                                            scenario.constraints.c_jac.', scenario.constraints.ceq_jac.', 'vars', {scenario.params.params, scenario.add_params.params});

	disp('Changing form of linear constraints to a more quick to evaluate form.')
	% It's actually quicker to simplify an expression, convert it to a function handle with matlabFunction, then evaluate the function
	% than to use subs (according to my benchmarks on scenario.linear_constraints.A with a moderately large number of shoots). The advantage
	% gets larger when we evaluate it multiple times (which often happens with fmincon-related functions).
	scenario.linear_constraints.A   = matlabFunction(simplify(scenario.linear_constraints.A),   'vars', {scenario.add_params.params});
	scenario.linear_constraints.b   = matlabFunction(simplify(scenario.linear_constraints.b),   'vars', {scenario.add_params.params});
	scenario.linear_constraints.Aeq = matlabFunction(simplify(scenario.linear_constraints.Aeq), 'vars', {scenario.add_params.params});
	scenario.linear_constraints.beq = matlabFunction(simplify(scenario.linear_constraints.beq), 'vars', {scenario.add_params.params});
end

% This runs a numerical optimization via fmincon using the scenario's objective and constraint functions
function scenario = traj_opt_fmincon(scenario)
	disp('Optimizing using fmincon')

	% We use these functions to call the actual functions; otherwise they will fail (as anonymous functions, they use deal to assign their outputs, which
	% fails when fmincon calls them with a different number of outputs)
	function [cost,gcost] = cost_fmincon(params)
		[cost,gcost] = scenario.fmincon.objective(params, scenario.add_params.value);
	end
	function [c,ceq,gc,gceq] = nonlcon_fmincon(params)
		[c,ceq,gc,gceq] = scenario.fmincon.nonlcon(params, scenario.add_params.value);
	end

	% Set up the optimization options
	disp('	Finalizing optimization options')
	default_options = optimset('Algorithm',  'sqp', ...
	                           'GradObj',    'on', ...
	                           'GradConstr', 'on');
	scenario.fmincon = declare_field(scenario.fmincon, 'options');
	scenario.fmincon.options = optimset(default_options, scenario.fmincon.options);

	% Set up the remainder of the necessary fields
	scenario                    = declare_field(scenario,                    'linear_constraints');
	scenario.params             = declare_field(scenario.params,             'lb');
	scenario.params             = declare_field(scenario.params,             'ub');

	% Evaluate the linear constraint functions with the given additional parameters
	disp('	Computing linear constraints with given additional parameters.')
	%A   = double(subs(scenario.linear_constraints.A,   scenario.add_params.params, scenario.add_params.value));
	%b   = double(subs(scenario.linear_constraints.b,   scenario.add_params.params, scenario.add_params.value));
	%Aeq = double(subs(scenario.linear_constraints.Aeq, scenario.add_params.params, scenario.add_params.value));
	%beq = double(subs(scenario.linear_constraints.beq, scenario.add_params.params, scenario.add_params.value));
	A   = scenario.linear_constraints.A(  scenario.add_params.value);
	b   = scenario.linear_constraints.b(  scenario.add_params.value);
	Aeq = scenario.linear_constraints.Aeq(scenario.add_params.value);
	beq = scenario.linear_constraints.beq(scenario.add_params.value);

	% Run the numerical optimization itself (3 times, since it doesn't always converge completely the first time).
	% In the future, this should be made arbitrary, and should terminate based on the return value, iter count, and other values from fmincon.
	disp('	Running fmincon...')
	scenario.minimum(:,1) = fmincon(@cost_fmincon, scenario.params.init_params, A, b, Aeq, beq, scenario.params.lb, scenario.params.ub, ...
		@nonlcon_fmincon, scenario.fmincon.options);
	scenario.minimum(:,2) = fmincon(@cost_fmincon, scenario.minimum(:,1),       A, b, Aeq, beq, scenario.params.lb, scenario.params.ub, ...
		@nonlcon_fmincon, scenario.fmincon.options);
	scenario.minimum(:,3) = fmincon(@cost_fmincon, scenario.minimum(:,2),       A, b, Aeq, beq, scenario.params.lb, scenario.params.ub, ...
		@nonlcon_fmincon, scenario.fmincon.options);
	disp('	fmincon-based optimization complete.')
end
