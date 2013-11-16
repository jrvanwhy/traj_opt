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
	if isfield(scenario, 'opt_fmincon')
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

	% Maximum iteration counts
	if ~isfield(scenario.opt_fmincon, 'max_seek_iters')
		scenario.opt_fmincon.max_seek_iters = 20;
	end
	if ~isfield(scenario.opt_fmincon, 'max_possible_iters')
		scenario.opt_fmincon.max_possible_iters = 5;
	end
	if ~isfield(scenario.opt_fmincon, 'max_found_iters')
		scenario.opt_fmincon.max_found_iters = 2;
	end

	% Reset the outputs
	scenario.opt_fmincon.soln           = [];
	scenario.opt_fmincon.fvals          = [];
	scenario.opt_fmincon.exit_flags     = [];
	scenario.opt_fmincon.output_structs = {};

	% Our initial point
	x0 = ones(scenario.param_maps.n_opt_params, 1);

	% Count of minimum-seeking iterations
	seek_iters = 0;

	% Count of "local minimum possible" iterations
	possible_iters = 0;

	% Count of "local minimum found" iterations
	found_iters = 0;

	% We loop here, calling fmincon repeatedly. Resetting the optimization seems to give better results; also, if rerun after a
	% "local minimum possible" output, it often finds a the local minimum to a higher accuracy. Also, this loop allows for more
	% robust behavior with respect to failing versus successful optimizations.
	done = false;
	while ~done
		% Call fmincon, retrieving variouus useful outputs
		[scenario.opt_fmincon.soln(:,end+1),       ...
		 scenario.opt_fmincon.fvals(end+1,1),      ...
		 scenario.opt_fmincon.exit_flags(end+1,1), ...
		 scenario.opt_fmincon.output_structs{end+1,1}] = ...
			fmincon(@obj_fcn, x0, [], [], [], [], [], [], @constr_fcn, options);

		% Switch based off the exit flag
		switch scenario.opt_fmincon.exit_flags(end)
			case -2
				% Not feasible!
				success = false;
				done    = true;
			case 0
				% Still seeking a minimum
				seek_iters = seek_iters + 1;

				% Check if we've iterated a lot
				if seek_iters > scenario.opt_fmincon.max_seek_iters
					success = false;
					done    = true;
				end
			case 1
				% Local minimum found!
				found_iters = found_iters + 1;

				% If we've hid our found iters limit, go ahead and return success
				if found_iters > scenario.opt_fmincon.max_found_iters
					success = true;
					done    = true;
				end
			case 2
				% Near a local minimum
				possible_iters = possible_iters + 1;

				% Check if we've stayed in this case for too long
				if possible_iters > scenario.opt_fmincon.max_possible_iters
					success = false;
					done    = true;
				end
		end

		% Update our initial value for the next iteration
		if ~done
			x0 = scenario.opt_fmincon.soln(:,end);
		end
	end

	% Copy over the solution as the final solution
	scenario.soln = scenario.opt_fmincon.soln(:,end);
end
