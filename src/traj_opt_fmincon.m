% This function runs a numerical optimization using MATLAB's fmincon function
% It runs the optimization for the given scenario, returning the scenario with optimization
% results as well as whether or not the optimization was successful
function [scenario,success] = traj_opt_fmincon(scenario, noopt_params, iterfcn)
	% These are the functions fmincon calls
	function [obj,gobj] = obj_fcn(params)
		% Iterate through the cost functions, summing them and their jacobians
		obj  = 0;
		jobj = sparse(1, scenario.n_params);
		for iterfcn = 1:numel(scenario.costs)
			cur_cost = scenario.costs{iterfcn};
			obj      = obj + sum(opt_eval_vecfcn(cur_cost, params));
			jobj     = jobj + sum(opt_eval_vecfcn_jac(cur_cost, params), 1);
		end

		% The "gradient" is actually the transpose of the jacobian
		gobj = jobj.';
	end
	function [c,ceq,gc,gceq] = constr_fcn(params)
		% Like for the cost, we need to iterate through the functions to call them all
		c   = [];
		jc  = sparse(0, scenario.n_params);
		for iterfcn = 1:numel(scenario.c)
			cur_c      = scenario.c{iterfcn};
			cur_c_val  = opt_eval_vecfcn(cur_c, params);
			c          = [c; cur_c_val(:)];
			cur_jc_val = opt_eval_vecfcn_jac(cur_c, params);
			jc         = [jc; cur_jc_val];
		end

		ceq  = [];
		jceq = sparse(0, scenario.n_params);
		for iterfcn = 1:numel(scenario.ceq)
			cur_ceq      = scenario.ceq{iterfcn};
			cur_ceq_val  = opt_eval_vecfcn(cur_ceq, params);
			ceq          = [ceq; cur_ceq_val(:)];
			cur_jceq_val = opt_eval_vecfcn_jac(cur_ceq, params);
			jceq         = [jceq; cur_jceq_val];
		end

		% The "gradient" is just the jacobian's transpose
		gc   = jc.';
		gceq = jceq.';
	end
	function hessian = hess_fcn(params, lambda)
		% Initialize hessian
		hessian = sparse(scenario.n_params, scenario.n_params);

		% Start by computing the hessian for each cost
		for itercost = 1:numel(scenario.costs)
			cur_cost = scenario.costs{itercost};
			hessian  = hessian + opt_eval_vecfcn_hess(cur_cost, params);
		end

		% Then add on the inequality constraint hessians
		lineqpos = 0;
		for iterc = 1:numel(scenario.c)
			cur_c         = scenario.c{iterc};
			cur_lambdas   = lambda.ineqnonlin(lineqpos+1:end);
			[c_hess,ladv] = opt_eval_vecfcn_hess(cur_c, params, cur_lambdas);
			hessian       = hessian + c_hess;
			lineqpos      = lineqpos + ladv;
		end

		% Last, add on the equality constraint hessians
		leqpos = 0;
		for iterceq = 1:numel(scenario.ceq)
			cur_ceq         = scenario.ceq{iterceq};
			cur_lambdas     = lambda.eqnonlin(leqpos+1:end);
			[ceq_hess,ladv] = opt_eval_vecfcn_hess(cur_ceq, params, cur_lambdas);
			hessian         = hessian + ceq_hess;
			leqpos          = leqpos + ladv;
		end
	end

	% Default iterfcn to empty
	if nargin < 3
		iterfcn = [];
	end

	% Optimization options
	options = optimset('Algorithm',       'interior-point', ...
	                   'Display',         'iter',           ...
	                   'Hessian',         'user-supplied',  ...
	                   'HessFcn',         @hess_fcn,        ...
	                   'GradObj',         'on',             ...
	                   'GradConstr',      'on',             ...
	                   'OutputFcn',       iterfcn);
	if isfield(scenario, 'fmincon') && isfield(scenario.fmincon, 'options')
		disp('	Using user-supplied options')
		options = optimset(options, scenario.fmincon.options);
	end
	scenario.opt_fmincon.options = options;

	% Maximum iteration counts
	if ~isfield(scenario.opt_fmincon, 'max_seek_iters')
		scenario.opt_fmincon.max_seek_iters = 100;
	end
	if ~isfield(scenario.opt_fmincon, 'max_possible_iters')
		scenario.opt_fmincon.max_possible_iters = 30;
	end
	if ~isfield(scenario.opt_fmincon, 'max_found_iters')
		scenario.opt_fmincon.max_found_iters = 0;
	end

	% Reset the outputs
	scenario.opt_fmincon.soln           = [];
	scenario.opt_fmincon.fvals          = [];
	scenario.opt_fmincon.exit_flags     = [];
	scenario.opt_fmincon.output_structs = {};

	% Our initial point
	x0 = scenario.x0;

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

				% If we've hit our found iters limit, go ahead and return success
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

% This function evalutes the hessian of the lagrangian of the given vecfcn.
% If you don't supply lambdas, it will default to all ones
function [hess_out,out_pos] = opt_eval_vecfcn_hess(fcn_in, params, lambdas)
	% Initialize lambdas if not given
	if nargin < 3
		lambdas = ones(fcn_in.n_vals * size(fcn_in.params, 2), 1);
	end

	% Initialize the hessian
	hess_out = sparse(numel(params), numel(params));

	% Set up initial output position
	out_pos = 0;

	% Iterate through each "equivalent function call", evaluating the hessian for each
	for iterout = 1:size(fcn_in.params, 2)
		% Parameters for this call
		call_params = fcn_in.params(:,iterout);

		% Grab the correct lambdas and update out_pos
		% Properly index into lambdas
		lambda_idxs     = out_pos+1:out_pos+fcn_in.n_vals;
		call_lambdas    = lambda_idxs.';
		call_lambdas(:) = lambdas(lambda_idxs);
		out_pos         = out_pos + fcn_in.n_vals;

		% Evaluate the hessian
		param_vals    = zeros(size(call_params));
		param_vals(:) = params(call_params);
		call_hess     = fcn_in.hess_fcn(param_vals, call_lambdas);

		% Update the sparse hessian
		hess_out(call_params, call_params) = hess_out(call_params, call_params) + call_hess;
	end
end

% This function evaluates the jacobian of the given vectorized
% function at a given point
function jval = opt_eval_vecfcn_jac(fcn_in, params)
	% Initialize the jacobian
	jval = sparse(size(fcn_in.params, 2) * fcn_in.n_vals, numel(params));

	% Iterate through each "equivalent function call",
	% evaluating the jacobian for each
	for iterout = 1:size(fcn_in.params, 2)
		% Parameters for this call
		call_params = fcn_in.params(:,iterout);

		% Evaluate the jacobian
		jac_params = zeros(size(call_params));
		jac_params(:) = params(call_params);
		call_jac = fcn_in.jac_fcn(jac_params);

		% Update the full jacobian values for this output
		rows = (size(call_jac, 1) * (iterout-1) + 1):(size(call_jac, 1) * iterout);
		cols = call_params;
		jval(rows,cols) = call_jac;
	end
end
