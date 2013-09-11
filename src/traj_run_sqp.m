% This function optimizes using a custom SQP implementation

function scenario = traj_run_sqp(scenario)
	scenario = traj_gen_sqp(scenario);
	scenario = traj_opt_sqp(scenario);
end

% Initialize all the values in scenario relevant to the SQP implementation
function scenario = traj_gen_sqp(scenario)
	% Check if we need to set up at all
	if isfield(scenario, 'sqp')
		scenario = scenario; % Make MATLAB happy
		return;
	end

	% We definitely need to generate the relevant functions
	disp('	Generating SQP structure.')

	% Some useful variables
	optvars    = scenario.params.params;
	add_params = scenario.add_params.params;
	obj_jac    = scenario.objective.jacobian;
	disp('		Calculating Hessian of objective.')
	obj_hess   = jacobian(obj_jac.', optvars);
	% The hessian must be symmetric... make it so
	disp('		Making Hessian of objective symmetric.')
	obj_hess   = (obj_hess + obj_hess.')/2;
	A          = scenario.linear_constraints.A;
	b          = scenario.linear_constraints.b;
	Aeq        = scenario.linear_constraints.Aeq;
	beq        = scenario.linear_constraints.beq;
	c          = scenario.constraints.c;
	c_jac      = scenario.constraints.c_jac;
	ceq        = scenario.constraints.ceq;
	ceq_jac    = scenario.constraints.ceq_jac;

	disp('		Generating hessian function.')
	scenario.sqp.obj_hess = matlabFunction(obj_hess,       'vars', {optvars, add_params});

	disp('		Generating objective jacobian function.')
	scenario.sqp.obj_jac  = matlabFunction(obj_jac,        'vars', {optvars, add_params});

	disp('		Generating A function.')
	% Apply corrections to re-center the linear constraints about the current point.
	scenario.sqp.A        = matlabFunction([A; c_jac],     'vars', {optvars, add_params});

	disp('		Generating b function.')
	scenario.sqp.b        = matlabFunction([b-A*optvars; -c],        'vars', {optvars, add_params});

	disp('		Generating Aeq function.')
	scenario.sqp.Aeq      = matlabFunction([Aeq; ceq_jac], 'vars', {optvars, add_params});

	disp('		Generating beq function.')
	scenario.sqp.beq      = matlabFunction([beq-Aeq*optvars; -ceq],    'vars', {optvars, add_params});

	disp('		Generating options.')
	scenario.sqp.options  = optimset('Algorithm', 'active-set', 'Display', 'off');
end

% Run the optimization!
function scenario = traj_opt_sqp(scenario)
	% Our current point
	x = scenario.params.init_params;

	% We use this for the diagnostic output
	iter = 0;

	% This is a filtered condition value -- it is used to correct the hessian
	% when the problem is found to be ill-conditioned.
	cond_filt = 1;

	% It's much easier to detect termination from within the loop.
	while true
		% Update iteration count
		iter = iter + 1;

		% Approximate the NLP problem with a QP problem.
		% Note that the origin in the approximation corresponds to the current
		% value of x in the real problem! Thus, the solution to the QP is a relative step.
		hess = scenario.sqp.obj_hess(x, scenario.add_params.value);
		grad = scenario.sqp.obj_jac( x, scenario.add_params.value);
		A    = scenario.sqp.A(       x, scenario.add_params.value);
		b    = scenario.sqp.b(       x, scenario.add_params.value);
		Aeq  = scenario.sqp.Aeq(     x, scenario.add_params.value);
		beq  = scenario.sqp.beq(     x, scenario.add_params.value);

		% This is a measure of how well-conditioned the problem is.
		condition = min(eig(hess));

		% Updated the filtered condition estimate
		cond_filt = cond_filt + (abs(condition) - cond_filt)/10;

		% Correct the hessian to be positive definite, if necessary
		if condition <= 0
			hess = hess + (cond_filt-condition) * eye(size(x, 1));
		end

		% Solve the QP subproblem to get the relative step to take
		% Note that we feed in the origin for x0 since the optimization is relative to x.
		% This is actually warm-starting quadprog!
		step = quadprog(hess, grad, A, b, Aeq, beq, [], [], zeros(size(x)), scenario.sqp.options);

		% Update x
		x    = x + step;

		% Diagnostics
		% Spit out a header every 40 iterations
		if mod(iter, 40) == 1
			disp('Iteration      Stepsize    Constraint Violation    Condition    Condition Estimate');
		end
		% And the diagnostic information itself
		fprintf('% 9u % 13f % 23f % 12f % 21f\n', iter, norm(step), max(abs([b(b < 0); beq])), condition, cond_filt);

		% If step was small enough, quit
		if dot(step, step) < 1e-20
			scenario.minimum = x;
			break;
		end

		% Debugging
		scenario.minimum = x;
		scenario = traj_gen_numerical(scenario);
		fplot(@(t) traj_get_state(scenario, t, 1), [0 scenario.num_duration])
		drawnow
	end
end
