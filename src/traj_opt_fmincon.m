% This function runs a numerical optimization using MATLAB's fmincon function

function [scenario,success] = traj_opt_fmincon(scenario, noopt_params)
	% Do the setup, if necessary
	scenario = setup_opt_fmincon(scenario);

	% Run the optimization
	[scenario, success] = run_opt_fmincon(scenario, noopt_params);
end

% This function does all the setup for an fmincon optimization run (if necessary)
function scenario = setup_opt_fmincon(scenario)
	% Exit if this has already been done
	if isfield(scenario, 'fmincon')
		scenario = scenario;
		return
	end

	% Diagnostic message for the user
	disp('	Doing setup for fmincon-based optimization.')

	% Create symbolic variables representing the optimization and non-optimized parameters
	opt_params   = sym('x', [scenario.param_maps.n_opt_params   1]);
	noopt_params = sym('n', [scenario.param_maps.n_noopt_params 1]);

	% Iterate through the cost functions to generate the objective
	objective = 0;
	for iter = 1:numel(scenario.optfcns.costs)
		objective = objective + scenario.optfcns.costs{iter}.fcn(opt_params, noopt_params);
	end
	% ... and simplify the objective
	objective = simplify(objective);

	% Iterate through the constraint functions to generate the constraint
	c   = sym([]);
	ceq = sym([]);
	for iter = 1:numel(scenario.optfcns.constraints)
		% Switch based on the type of constraint
		if scenario.optfcns.constraints{iter}.con_type == '<='
			c   = [c;   reshape(scenario.optfcns.constraints{iter}.fcn(opt_params, noopt_params), [], 1)];
		else
			ceq = [ceq; reshape(scenario.optfcns.constraints{iter}.fcn(opt_params, noopt_params), [], 1)];
		end
	end
	% now simplify the constraints
	c   = simplify(c);
	ceq = simplify(ceq);

	% For now, we'll just take the gradients of the full optimization and constraint functions.
	% Later, we'll incorporate them earlier in the line, but that'll take a lot more work, and has not yet
	% been shown to be necessary. Besides, if I do MPC stuff, fmincon-based optimization may become a thing of the past...
	scenario.opt_fmincon.obj_fcn     = matlabFunction(objective, 'vars', {opt_params, noopt_params});
	scenario.opt_fmincon.c_fcn       = matlabFunction(c,         'vars', {opt_params, noopt_params});
	scenario.opt_fmincon.ceq_fcn     = matlabFunction(ceq,       'vars', {opt_params, noopt_params});
	scenario.opt_fmincon.gradobj_fcn = matlabFunction(simplify(jacobian(objective, opt_params).'), 'vars', {opt_params, noopt_params});
	scenario.opt_fmincon.gradc_fcn   = matlabFunction(simplify(jacobian(c,         opt_params).'), 'vars', {opt_params, noopt_params});
	scenario.opt_fmincon.gradceq_fcn = matlabFunction(simplify(jacobian(ceq,       opt_params).'), 'vars', {opt_params, noopt_params});

	% Clean up symbolic variables
	disp('	Cleaning up symbolic variables')
	syms x n clear
end

% This function runs the optimization for the given scenario, returning the scenario with optimization
% results in at as well as whether or not the optimization was successful
function [scenario,success] = run_opt_fmincon(scenario, noopt_params)
	% We use proxy functions both to deal with noopt_params and to
	% deal with the varying number of output arguments fmincon expects.
	function [obj,gobj] = obj_fcn(params)
		obj = scenario.opt_fmincon.obj_fcn(params, noopt_params);
		if nargout > 1
			gobj = scenario.opt_fmincon.gradobj_fcn(params, noopt_params);
		end
	end
	function [c,ceq,gc,gceq] = constr_fcn(params)
		c = scenario.opt_fmincon.c_fcn(params, noopt_params);
		if nargout > 1
			ceq = scenario.opt_fmincon.ceq_fcn(params, noopt_params);
			if nargout > 2
				gc = scenario.opt_fmincon.gradc_fcn(params, noopt_params);
				if nargout > 3
					gceq = scenario.opt_fmincon.gradceq_fcn(params, noopt_params);
				end
			end
		end
	end

	% Optimization options
	options = optimset('Algorithm', 'interior-point', ...
	                   'Display',   'iter',           ...
	                   'GradObj',   'on',             ...
	                   'GradConstr', 'on');
	if isfield(scenario.opt_fmincon, 'options')
		options = optimset(scenario.opt_fmincon.options, options);
	end
	scenario.opt_fmincon.options = options;

	% Call fmincon
	scenario.opt_fmincon.soln = fmincon(@obj_fcn, ones(scenario.param_maps.n_opt_params, 1), [], [], [], [], [], [], @constr_fcn, options);

	% Copy over the solution as the final solution
	scenario.soln = scenario.opt_fmincon.soln;

	% For now, simply return success
	success = true;
end
