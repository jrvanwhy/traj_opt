% This function runs a numerical optimization using MATLAB's fmincon function

function [scenario,success] = traj_opt_fmincon(scenario, noopt_params)
	% Run the optimization
	[scenario, success] = run_opt_fmincon(scenario, noopt_params);
end

% This function runs the optimization for the given scenario, returning the scenario with optimization
% results in at as well as whether or not the optimization was successful
function [scenario,success] = run_opt_fmincon(scenario, noopt_params)
	% We use proxy functions both to deal with noopt_params and to
	% deal with the varying number of output arguments fmincon expects.
	function [obj,gobj] = obj_fcn(params)
		obj  = scenario.optfcns.cost(params);
		% The "gradient" is actually the transpose of the jacobian
		gobj = scenario.optfcns.jcost(params).';
	end
	function [c,ceq,gc,gceq] = constr_fcn(params)
		c    = scenario.optfcns.c(params);
		ceq  = scenario.optfcns.ceq(params);
		% Again, fmincon's "gradients" are just transposes of the jacobian
		gc   = sparse(scenario.optfcns.jc.j,           ...
		              scenario.optfcns.jc.i,           ...
		              scenario.optfcns.jc.s(params),   ...
		              scenario.optfcns.jc.n,           ...
		              scenario.optfcns.jc.m);
		gceq = sparse(scenario.optfcns.jceq.j,         ...
		              scenario.optfcns.jceq.i,         ...
		              scenario.optfcns.jceq.s(params), ...
		              scenario.optfcns.jceq.n,         ...
		              scenario.optfcns.jceq.m);
	end
	function hessian = hess_fcn(params, lambda)
		hessian = sparse(scenario.optfcns.hess.i,                                             ...
		                 scenario.optfcns.hess.j,                                             ...
		                 scenario.optfcns.hess.s(params, lambda.ineqnonlin, lambda.eqnonlin), ...
		                 scenario.optfcns.hess.mn,                                            ...
		                 scenario.optfcns.hess.mn);
	end

	% Optimization options
	options = optimset('Algorithm',       'interior-point', ...
	                   'Display',         'iter',           ...
	                   'GradObj',         'on',             ...
	                   'GradConstr',      'on');
	% Use a user-specified Hessian if available. Otherwise, bfgs is probably better
	if isfield(scenario.optfcns, 'hess')
		disp('	Using sparse user-specified Hessian')
		options = optimset(options, 'Hessian', 'user-supplied', 'HessFcn', @hess_fcn);
	end
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
		scenario.opt_fmincon.max_found_iters = 2;
	end

	% Reset the outputs
	scenario.opt_fmincon.soln           = [];
	scenario.opt_fmincon.fvals          = [];
	scenario.opt_fmincon.exit_flags     = [];
	scenario.opt_fmincon.output_structs = {};

	% Our initial point
	x0 = ones(scenario.n_params, 1);

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
