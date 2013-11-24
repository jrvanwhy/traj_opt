% This is a benchmarking function for the optimizer. It uses the
% slip hop example, but that example was encoded into this file so each part
% may be tested separately.

function results = benchmark(cutoff_time, ode_type)
	% The current number of intervals
	cur_intrvls = 2; % Start at 2, as that's the lowest value for which the optimization converges.

	% Initialize the results
	results.t_phase    = [];
	results.t_scenario = [];
	results.t_opts     = [];
	results.phases     = {};
	results.scenarios  = {};
	results.opts       = {};
	results.fvals      = [];

	% Iterate through the intervals until we reach the cutoff time
	while sum(results.t_phase) + sum(results.t_scenario) + sum(results.t_opts) < cutoff_time
		% Do phase setup
		tic;
		results.phases{1,end+1} = traj_setup_dynamic_phase('stance',                          ...
						                   @dynamic_system,                   ...
						                   {'position', 'velocity', 'force'}, ...
						                   {'actuator_rate'},                 ...
						                   cur_intrvls,                       ...
						                   [],                                ...
						                   [],                                ...
						                   ode_type);
		results.t_phase(1,end+1) = toc;

		% Do scenario setup
		tic;
		results.scenarios{1,end+1} = traj_setup_scenario(@scenario_fcn, results.phases{end});
		results.t_scenario(1,end+1) = toc;

		% Do the optimization
		tic;
		results.opts{1,end+1}   = traj_optimize(results.scenarios{end});
		results.t_opts(1,end+1) = toc;
		results.fvals(1,end+1)  = results.opts{1,end}.opt_fmincon.fvals(end);

		% Increment current interval count
		cur_intrvls = cur_intrvls + 1;
	end
end

%% Let's return the completed scenario structure, so the user may inspect the results
%% Also, the number of "intervals" for the dynamics representation (equivalent to timesteps for an ODE solver)
%% is a parameter of this function.
%function scenario = slip_hop(n_intervals)
%	% Here, we create a dynamic phase representing the stance phase. The first parameter is the name of the phase (optional), the second
%	% should be the "dynamic system" function, which is of the form [dx,...] = dynsys(x, u, add_params, noopt_params)
%	% and defines costs and constraints throughout this dynamic phase. The third parameter defines the state
%	% for this phase, and the fourth defines the input. The fifth parameter is the number of time intervals in the
%	% discrete approximation to the trajectory and input. The sixth and seventh are optional, and specify additional
%	% optimization parameters and non-optimized parameters, neither of which we need for this scenario.
%	phase = traj_setup_dynamic_phase('stance',                          ...
%	                                 @dynamic_system,                   ...
%	                                 {'position', 'velocity', 'force'}, ...
%	                                 {'actuator_rate'},                 ...
%	                                 n_intervals);
%
%	% This sets up the scenario. We pass a name (optional, not used here),  the scenario function
%	% (optional), and each phase for the scenario (in order).
%	% If multiple phases with the same name are passed, this will automatically append "_N" to their name, where
%	% N is an increasing positive integer.
%	scenario = traj_setup_scenario(@scenario_fcn, phase);
%
%	% Run the optimization
%	% This returns the scenario and a success flag (also stored in the scenario)
%	% that is true if it was successful and false otherwise.
%	[scenario,success] = traj_optimize(scenario);
%
%	% Check if it failed -- if so, run the constraint feasibility analyzer
%	if ~success
%		% It failed -- run the analyzer
%		traj_analyze_feasibility(scenario);
%	end
%end

% This is the dynamic system function. It does most of the work of specifying the dynamic behavior
% of the inverted pendulum and constraints on its dynamic motion.
% Since add_params and noopt_params is not necessary for this dynamic system, we do not need to make
% them parameters.
function [dx,cost,powerfcn] = dynamic_system(x, u)
	% Some basic system parameters
	mass         = 59.9; % Mass of the SLIP in kg
	spring_const = 6543; % Spring constant on N/m
	g            = 9.81; % Gravitational acceleration, m/s²

	% Let's give the relevant variables better names
	position = x(1);
	velocity = x(2);
	force    = x(3);
	act_rate = u(1);

	% The derivative of position is velocity; the derivative of velocity is calculated from force
	% and gravitational acceleration. The derivative of force is calculated from the spring constant
	% and actuator rate.
	dx = [ velocity
	       force/mass - g
	       spring_const * (act_rate - velocity)];

	% This is power, used as for our cost function
	power    = act_rate * force;
	powerfcn = traj_create_function('power', power);

	% Set up our cost (unsigned work)
	c    = 1;
	cost = traj_create_cost('work', sqrt(power^2 + c^2)-c);
end

% This is the scenario function. It takes in a scenario that is populated with symbolic values, and imposes additional constraints
% and costs upon the scenario.
function [start_state_con,end_state_con] = scenario_fcn(scenario)
	% Touchdown parameters
	g         = 9.81;                % Gravitational accel, m/s²
	td_hgt    = .8;                  % Height at touchdown
	drop_dist = .2;                  % Drop height, m
	td_spd    = sqrt(2*g*drop_dist); % Speed at touchdown

	% Constrain the starting state to match our given touchdown conditions, and our ending state to be its "opposite".
	start_state_con = traj_create_constraint('Starting state', scenario.phases{1}.states(:,1),   '=', [td_hgt; -td_spd; 0]);
	end_state_con   = traj_create_constraint('Ending state',   scenario.phases{1}.states(:,end), '=', [td_hgt;  td_spd; 0]);
end
