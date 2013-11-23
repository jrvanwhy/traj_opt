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

	% Generate the final objecttive and constraint functions
	scenario = gen_optfcns(scenario, scenario_fcn);
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

% This does the remainder of the necessary setup for optimization
function scenario = gen_optfcns(scenario, scenario_fcn)
	% Set parameter count
	scenario.n_params = scenario.param_map(2,end);

	disp('	Creating symbolic optimization parameters')
	opt_params = sym('x', [scenario.n_params 1]);

	% Add on remaining costs and constraints
	scenario = process_scenario_fcn(scenario, scenario_fcn, opt_params);

	disp('	Creating final cost expression')
	cost_expr = sym(0);
	% Iterate through the phases to sum costs
	for iter_phase = 1:numel(scenario.phases)
		phase = scenario.phases{iter_phase};
		disp(['		Adding costs from phase ''' phase.names.phase ''''])
		phase_params = opt_params(scenario.param_map(1,iter_phase):scenario.param_map(2,iter_phase));

		cost_expr = cost_expr + phase.cost(phase_params);
	end
	% Iterate through the scenario function costs
	for iter_scen_cost = 1:numel(scenario.scenfun.costs)
		cost = scenario.scenfun.costs{iter_scen_cost};
		disp(['		Adding cost ''' cost.name ''''])

		cost_expr = cost_expr + cost.fcn(opt_params);
	end
	disp('	Simplifying final cost')
	cost_expr = simplify(cost_expr);
	disp('	Converting final cost to anonymous function')
	scenario.optfcns.cost = matlabFunction(cost_expr, 'vars', {opt_params});
	disp('	Creating cost gradient expression')
	gcost_expr = jacobian(cost_expr, opt_params).';
	disp('	Converting cost gradient to anonymous function')
	scenario.optfcns.gcost = matlabFunction(gcost_expr, 'vars', {opt_params});

	disp('	Creating final constraint expressions')
	c_expr                  = sym([]);
	ceq_expr                = sym([]);
	gc_s_expr               = sym([]);
	gceq_s_expr             = sym([]);
	scenario.optfcns.gc.i   = [];
	scenario.optfcns.gc.j   = [];
	scenario.optfcns.gc.m   = 0;
	scenario.optfcns.gc.n   = 0;
	scenario.optfcns.gceq.i = [];
	scenario.optfcns.gceq.j = [];
	scenario.optfcns.gceq.m = 0;
	scenario.optfcns.gceq.n = 0;
	% Iterate through the phases, processing their constraints
	for iter_phase = 1:numel(scenario.phases)
		phase = scenario.phases{iter_phase};
		disp(['		Adding constraints from phase ''' phase.names.phase ''''])
		phase_params = opt_params(scenario.param_map(1,iter_phase):scenario.param_map(2,iter_phase));

		c_expr                  = [c_expr;                   phase.c(phase_params)               ];
		ceq_expr                = [ceq_expr;                 phase.ceq(phase_params)             ];
		gc_s_expr               = [gc_s_expr;                phase.gc.s(phase_params)            ];
		gceq_s_expr             = [gceq_s_expr;              phase.gceq.s(phase_params)          ];
		scenario.optfcns.gc.i   = [scenario.optfcns.gc.i;    phase.gc.i+scenario.optfcns.gc.m    ];
		scenario.optfcns.gc.j   = [scenario.optfcns.gc.j;    phase.gc.j+scenario.optfcns.gc.n    ];
		scenario.optfcns.gc.m   =  scenario.optfcns.gc.m   + phase.gc.m;
		scenario.optfcns.gc.n   =  scenario.optfcns.gc.n   + phase.gc.n;
		scenario.optfcns.gceq.i = [scenario.optfcns.gceq.i;  phase.gceq.i+scenario.optfcns.gceq.m];
		scenario.optfcns.gceq.j = [scenario.optfcns.gceq.j;  phase.gceq.j+scenario.optfcns.gceq.n];
		scenario.optfcns.gceq.m =  scenario.optfcns.gceq.m + phase.gceq.m;
		scenario.optfcns.gceq.n =  scenario.optfcns.gceq.n + phase.gceq.n;
	end
	% Iterate through scenario function constraints
	for iter_scen_con = 1:numel(scenario.scenfun.constraints)
		con = scenario.scenfun.constraints{iter_scen_con};
		disp(['		Adding constraint ''' con.name ''''])

		con_expr                    = con.fcn(opt_params);
		gcon_expr                   = jacobian(con_expr, opt_params).';
		gcon_n                      = size(gcon_expr, 2);
		[gcon_i,gcon_j,gcon_s_expr] = find(gcon_expr);

		if con.con_type == '<='
			c_expr                  = [c_expr;                   con_expr                      ];
			gc_s_expr               = [gc_s_expr;                gcon_s_expr                   ];
			scenario.optfcns.gc.i   = [scenario.optfcns.gc.i;    gcon_i                        ];
			scenario.optfcns.gc.j   = [scenario.optfcns.gc.j;    gcon_j+scenario.optfcns.gc.n  ];
			scenario.optfcns.gc.n   =  scenario.optfcns.gc.n   + gcon_n;
		else
			ceq_expr                = [ceq_expr;                 con_expr                      ];
			gceq_s_expr             = [gceq_s_expr;              gcon_s_expr                   ];
			scenario.optfcns.gceq.i = [scenario.optfcns.gceq.i;  gcon_i                        ];
			scenario.optfcns.gceq.j = [scenario.optfcns.gceq.j;  gcon_j+scenario.optfcns.gceq.n];
			scenario.optfcns.gceq.n =  scenario.optfcns.gceq.n + gcon_n;
		end
	end
	disp('	Simplifyng inequality constraints')
	c_expr   = simplify(c_expr);

	disp('	Simplifyng equality constraints')
	ceq_expr = simplify(ceq_expr);

	disp('	Simplifying inequality constraint gradient')
	gc_s_expr   = simplify(gc_s_expr);

	disp('	Simplifying equality constraint gradient')
	gceq_s_expr = simplify(gceq_s_expr);

	disp('	Creating inequality constraint anonymous function')
	scenario.optfcns.c   = matlabFunction(c_expr,   'vars', {opt_params});

	disp('	Creating equality constraint anonymous function')
	scenario.optfcns.ceq = matlabFunction(ceq_expr, 'vars', {opt_params});

	disp('	Creating inequality constraint gradient anonymous function')
	scenario.optfcns.gc.s   = matlabFunction(gc_s_expr,   'vars', {opt_params});

	disp('	Creating equality constraint gradient anonymous function')
	scenario.optfcns.gceq.s = matlabFunction(gceq_s_expr, 'vars', {opt_params});

	disp('	Cleaning up symbolic variables')
	syms opt_params clear
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
end

% Handles everything necessary for the scenario function
function scenario = process_scenario_fcn(scenario, scenario_fcn, opt_params)
	% Check if a scenario function was given. If not, return scenario.
	if isempty(scenario_fcn)
		scenario = scenario;
		return
	end

	disp('	Setting up to call scenario function')

	% Add the symbolic outputs for each function to the scenario
	scenario = traj_eval_funcs(scenario, opt_params, []);

	% Call the scenario function, capturing all output arguments
	disp('	Calling scenario function')
	[scenario_varouts{1:nargout(scenario_fcn)}] = scenario_fcn(scenario);

	% Initialize scenario function fields
	scenario.scenfun.constraints = {};
	scenario.scenfun.costs       = {};

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

				% Convert the constraint function to an anonymous function, then append it to the
				% constraints field.
				output.fcn = matlabFunction(simplify(output.fcn), 'vars', ...
					{opt_params});
				scenario.scenfun.constraints{end+1} = output;

			case 'cost'
				% Progress message so the user doesn't think it's frozen
				disp(['		Processing cost ''' output.name ''''])

				% Convert the cost function to an anonymous function with appropriate inputs, then
				% append it to the costs field
				output.fcn = matlabFunction(output.fcn, 'vars', ...
					{opt_params});
				scenario.scenfun.costs{end+1} = output;

			otherwise
				% Spit out an error -- this is not an acceptable structure type
				error(['Structure type ''' output.struct_type ''' not a valid return type for the dynamic system function.'])
		end
	end

	disp(['	Cleaning up symbolic variables'])
end
