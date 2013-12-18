% This function runs a numerical optimization using MATLAB's fmincon function
% It runs the optimization for the given scenario, returning the scenario with optimization
% results as well as whether or not the optimization was successful
function [scenario,success] = traj_opt_fmincon(scenario, noopt_params, iterfcn)
	% These are the functions fmincon calls
	function [obj,gobj] = obj_fcn(params)
		% Iterate through the cost functions, summing them and generating their jacobian
		obj    = 0;
		jobj_j = [];
		jobj_s = [];
		for iterfcn = 1:numel(scenario.costs)
			cur_cost = scenario.costs{iterfcn};
			obj      = obj + sum(opt_eval_vecfcn(cur_cost, params));
			[~,cur_jobj_j,cur_jobj_s] = opt_eval_vecfcn_jac(cur_cost, params);
			jobj_j = [jobj_j; cur_jobj_j];
			jobj_s = [jobj_s; cur_jobj_s];
		end

		gobj = sparse(jobj_j, ones(numel(jobj_j), 1), jobj_s, numel(params), 1);
	end
	function [c,ceq,gc,gceq] = constr_fcn(params)
		% Like for the cost, we need to iterate through the functions to call them all
		c    = [];
		jc_i = [];
		jc_j = [];
		jc_s = [];
		for iterfcn = 1:numel(scenario.c)
			cur_c      = scenario.c{iterfcn};
			cur_c_val  = opt_eval_vecfcn(cur_c, params);
			c          = [c; cur_c_val(:)];
			[cur_jc_i,cur_jc_j,cur_jc_s] = opt_eval_vecfcn_jac(cur_c, params);
			jc_i = [jc_i; cur_jc_i];
			jc_j = [jc_j; cur_jc_j];
			jc_s = [jc_s; cur_jc_s];
		end

		ceq    = [];
		jceq_i = [];
		jceq_j = [];
		jceq_s = [];
		for iterfcn = 1:numel(scenario.ceq)
			cur_ceq      = scenario.ceq{iterfcn};
			cur_ceq_val  = opt_eval_vecfcn(cur_ceq, params);
			ceq          = [ceq; cur_ceq_val(:)];
			[cur_jceq_i,cur_jceq_j,cur_jceq_s] = opt_eval_vecfcn_jac(cur_ceq, params);
			jceq_i = [jceq_i; cur_jceq_i];
			jceq_j = [jceq_j; cur_jceq_j];
			jceq_s = [jceq_s; cur_jceq_s];
		end

		% The "gradient" is just the jacobian's transpose
		gc   = sparse(jc_j,   jc_i,   jc_s,   numel(params), numel(c));
		gceq = sparse(jceq_j, jceq_i, jceq_s, numel(params), numel(ceq));
	end
	function hessian = hess_fcn(params, lambda)
		% Initialize hessian arrays
		hess_i = [];
		hess_j = [];
		hess_s = [];

		% Start by computing the hessian for each cost
		for itercost = 1:numel(scenario.costs)
			cur_cost = scenario.costs{itercost};
			[cur_hess_i,cur_hess_j,cur_hess_s] = opt_eval_vecfcn_hess(cur_cost, params);
			hess_i = [hess_i; cur_hess_i];
			hess_j = [hess_j; cur_hess_j];
			hess_s = [hess_s; cur_hess_s];
		end

		% Then add on the inequality constraint hessians
		lineqpos = 0;
		for iterc = 1:numel(scenario.c)
			cur_c         = scenario.c{iterc};
			cur_lambdas   = lambda.ineqnonlin(lineqpos+1:end);
			[cur_hess_i,cur_hess_j,cur_hess_s] = opt_eval_vecfcn_hess(cur_c, params, cur_lambdas);
			hess_i = [hess_i; cur_hess_i];
			hess_j = [hess_j; cur_hess_j];
			hess_s = [hess_s; cur_hess_s];
			lineqpos      = lineqpos + numel(cur_c.fcn.const_outs);
		end

		% Last, add on the equality constraint hessians
		leqpos = 0;
		for iterceq = 1:numel(scenario.ceq)
			cur_ceq         = scenario.ceq{iterceq};
			cur_lambdas     = lambda.eqnonlin(leqpos+1:end);
			[cur_hess_i,cur_hess_j,cur_hess_s] = opt_eval_vecfcn_hess(cur_ceq, params, cur_lambdas);
			hess_i = [hess_i; cur_hess_i];
			hess_j = [hess_j; cur_hess_j];
			hess_s = [hess_s; cur_hess_s];
			leqpos          = leqpos + numel(cur_ceq.fcn.const_outs);
		end

		% Generate the hessian
		hessian = sparse(hess_i, hess_j, hess_s, numel(params), numel(params));
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
	                   'DerivativeCheck', 'off',            ...
	                   'FinDiffType',     'central',        ...
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
		scenario.opt_fmincon.max_found_iters = 1;
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

% This function evaluates the hessian of the lagrangian of the given vecfcn.
% If you don't supply lambdas, it will default to all ones
% It returns the i, j, and s values for the matrix (it's sparse)
function [hess_i,hess_j,hess_s] = opt_eval_vecfcn_hess(fcn_in, params, lambdas)
	% Initialize lambdas if not given
	if nargin < 3
		lambdas = ones(numel(fcn_in.fcn.const_outs), 1);
	end

	% Grab the right parameters in the right shape
	vecparams    = zeros(size(fcn_in.param_nums));
	vecparams(:) = params(fcn_in.param_nums);

	% Copy over the i and j vectors
	hess_i = fcn_in.hess_i;
	hess_j = fcn_in.hess_j;

	% Initialize hess_s using the constant outputs matrix
	hess_s = fcn_in.hess_s.const_outs;

	% Reshape the lambdas to be the right shape
	lambdas_new    = zeros(size(fcn_in.fcn.const_outs));
	lambdas_new(:) = lambdas(1:numel(lambdas_new));

	% Add in the nonconstant s values
	if nnz(fcn_in.hess_s.nc_outs) > 0
		hess_s(fcn_in.hess_s.nc_outs,:) = fcn_in.hess_s.nc_fcn(vecparams, lambdas_new);
	end

	% Reshape s into a column vector
	hess_s = hess_s(:);
end

% This function evaluates the jacobian of the given vectorized
% function at a given point
% It returns the i, j, and s values for the matrix
function [jac_i,jac_j,jac_s] = opt_eval_vecfcn_jac(fcn_in, params)

	% Grab the right parameters in the right shape
	vecparams    = zeros(size(fcn_in.param_nums));
	vecparams(:) = params(fcn_in.param_nums);

	% Copy over the i and j vectors
	jac_i = fcn_in.jac_i;
	jac_j = fcn_in.jac_j;

	% Initialize jac_s using the constant outputs matrix
	jac_s = fcn_in.jac_s.const_outs;

	% Add in the nonconstant s values
	if nnz(fcn_in.jac_s.nc_outs) > 0
		jac_s(fcn_in.jac_s.nc_outs,:) = fcn_in.jac_s.nc_fcn(vecparams);
	end

	% Reshape s into a column vector
	jac_s = jac_s(:);
end
