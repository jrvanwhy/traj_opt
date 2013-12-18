% This is sample code to demonstrate usage of additional optimization parameters.
% It solves a simple optimization problem with a known exact solution.
%
% The problem setup is as such: There is a car. At time 0, it is stationary at the origin.
% A time-varying force may be applied to this car. After some amount of time, the car must be stationary
% at position 1. This time should be minimized without exceeding the minimum and maximum force limitations.
%
% The input here (u) consists only of the force applied on the car. The state is a member of R²,
% and is of the following form:
% x = [ position ]
%     [ velocity ]
% For testing purposes, the input is represented as a polynomial
%
% The only scenario-imposed constraints are the force limits, starting state, and ending state.

% Let's return the completed scenario structure, so the user may inspect the results
% Also, the number of "intervals" for the dynamics representation (equivalent to timesteps for an ODE solver)
% is a parameter of this function.
function scenario = car_poly(n_intervals, poly_deg)
	% Define the additional parameters -- these are the fourier coefficients
	% (and a duration parameter for scaling, since time-variant systems aren't supported yet)
	add_params = {};
	for iter = 0:poly_deg
		add_params{end+1} = ['poly' num2str(iter)];
	end

	% Here, we create a dynamic phase representing the car's trip. The first parameter is the name of the phase (optional), the second
	% should be the "dynamic system" function, which is of the form [dx,...] = dynsys(x, u, add_params, noopt_params)
	% and defines costs and constraints throughout this dynamic phase. The third parameter defines the state
	% for this phase, and the fourth defines the input. The fifth parameter is the number of time intervals in the
	% discrete approximation to the trajectory and input. The sixth and seventh are optional, and specify additional
	% optimization parameters and non-optimized parameters, neither of which we need for this scenario.
	% The eighth parameter, if given, specifies the type of ODE solver to use.
	phase = traj_setup_dynamic_phase('transit',                        ...
	                                 @dynamic_system,                  ...
	                                 {'position', 'velocity', 'time'}, ...
	                                 [],                               ...
	                                 n_intervals,                      ...
	                                 add_params,                       ...
	                                 [],                               ...
	                                 'midpoint');

	% This sets up the scenario. We pass a name (optional, not used here),  the scenario function
	% (optional), and each phase for the scenario (in order).
	% If multiple phases with the same name are passed, this will automatically append "_N" to their name, where
	% N is an increasing positive integer.
	scenario = traj_setup_scenario(@scenario_fcn, phase);

	% Run the optimization
	% This returns the scenario and a success flag (also stored in the scenario)
	% that is true if it was successful and false otherwise.
	clf
	[scenario,success] = traj_optimize(scenario, [], @plot_input);

	% Check if it failed -- if so, run the constraint feasibility analyzer
	if ~success
		% It failed -- run the analyzer
		traj_analyze_feasibility(scenario);
	end
end

% This is the dynamic system function. It does most of the work of specifying the dynamic behavior
% of the car and constraints on its dynamic motion.
% Since add_params and noopt_params is not necessary for this dynamic system, we do not need to make
% them parameters.
function [dx,min_force_con,max_force_con,inputfcn] = dynamic_system(x, a)
	% Set the input as a function of the additional parameters and time (fourier series)
	time    = x(3);
	polyvar = sym(1);
	u       = sym(0);
	for iter = 1:numel(a)
		u       = u + a(iter) * polyvar / factorial(iter);
		polyvar = polyvar * time;
	end

	% The derivative of position is velocity; the derivative of velocity is the input, since
	% the mass is a nondimensional 1.
	dx = [ x(2)
	       u
	       1 ];

	% Make sure the absolute value of the force does not exceed 1
	min_force_con = traj_create_constraint('Min Force', u, '>=', -1);
	max_force_con = traj_create_constraint('Max Force', u, '<=',  1);

	% Create a function for the input because it's not so obvious anymore
	inputfcn = traj_create_function('force', u);
end

% This is the scenario function. It takes in a scenario that is populated with symbolic values, and imposes additional constraints
% and costs upon the scenario.
function [cost,start_state_con,end_state_con] = scenario_fcn(scenario)
	% Generate a cost (the only cost) proportional to the duration of the trajectory.
	cost = traj_create_cost('Duration', scenario.phases{1}.duration);

	% Constrain the starting state to be at the origin and the ending state to be
	% a stationary car at position 1
	start_state_con = traj_create_constraint('Starting state', scenario.phases{1}.states(:,1),     '=', [0; 0; 0]);
	end_state_con   = traj_create_constraint('Ending state',   scenario.phases{1}.states(1:2,end), '=', [1; 0]);
end

% This function simply plots the input
function plot_input(scenario)
	plot(scenario.phases{1}.t_force, scenario.phases{1}.force)
	drawnow
end