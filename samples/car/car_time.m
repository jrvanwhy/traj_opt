% This is sample code to demonstrate usage of the optimization framework.
% It also serves as example code in developing the new API.
% It solves a simple optimization problem with a known exact solution.
%
% The problem setup is as such: There is a car. At time 0, it is stationary at the origin.
% A time-varying force may be applied to this car. After some amount of time, the cart must be stationary
% at position 1. This time should be minimized without exceeding the maximum and minimum force limitations.
%
% The input here (u) consists only of the force applied on the car. The state is a member of R²,
% and is of the following form:
% x = [ position ]
%     [ velocity ]
%
% The only constraints are the force limits, starting state, and ending state.

% Let's return the completed scenario structure, so the user may inspect the results
% Also, the number of "intervals" for the dynamics representation (equivalent to timesteps for an ODE solver)
% is a parameter of this function.
function scenario = car_time(n_intervals)
	% This creates a state structure. We specify names for each element in the state.
	% The state will be a column vector, even if the names array isn't. The
	% number of elements in the names array determines the size of the state.
	% If a string is passed in, the state will be of size 1
	state = traj_create_state('position', 'velocity');

	% This is a similar to traj_create_state(), but for input instead.
	input = traj_create_input('force');

	% This is our dynamics function. It takes in the state and input and returns the derivative
	% of the state. It gets called by the optimization framework (during the call to traj_create_dynamics())
	function dx = car_dynamics(x, u)
		% We'll use a linear representation of the system here
		A = [ 0 1
		      0 0 ];

		B = [ 0
		      1 ];

		% Specify the derivative of the state
		dx = A * x + B * u;
	end

	% Create the dynamics structure
	dynamics = traj_create_dynamics(@car_dynamics, state, input);

	% Here, we create a dynamic phase. We must specify the dynamics
	% in use and may specify the number of intervals for this phase.
	% Additionally, we may create a "cost" function via traj_create_cost()
	% or constraints and add them to the resulting dynamic phase.
	phase = traj_create_phase(dynamics, n_intervals);

	% Now we need to constrain the starting and ending states to be equal to the origin
	% and one unit in the positive direction, respectively. To do this, we first need to be able
	% to access the states. We can get a vector of all the states using traj_get_states(phase)
	states = traj_get_states(phase);

	% Create the constraint for the starting state. Set it to the origin. traj_create_constraint()
	% takes three parameters. The 2nd parameter sets the constraint type ('='/'==', '<=', or '>='), and the other
	% two are the "sides" of the equation or inequality.
	initial_state_constraint = traj_create_constraint(states(:,1), '=', [ 0
	                                                                      0 ]);

	% Name this constraint. traj_add_name() can take in essentially any optimizer structure and
	% give it a name.
	initial_state_constraint = traj_add_name(initial_state_constraint, 'Initial state');

	% Similarly, create an inequality constraint for the ending state and name it
	ending_state_constraint = traj_create_constraint(states(:,end), '=', [ 1
	                                                                       0 ]);
	ending_state_constraint = traj_add_name(ending_state_constraint, 'Ending state');

	% Add the constraints to the phase
	% traj_add_constraint() adds the constraint to the given phase or scenario
	phase = traj_add_constraint(phase, initial_state_constraint);
	phase = traj_add_constraint(phase, ending_state_constraint);

	% Create the scenario structure itself. Go ahead and add the phase to the scenario during initialization.
	scenario = traj_create_scenario(phase);

	% Retrieve the duration of the scenario. In this case, traj_get_duration()
	% returns the sum of the durations of the phases. traj_get_duration() may also be used
	% on an individual phase to retrieve just its duration.
	duration = traj_get_duration(scenario);

	% Add on the cost. In our case, we're trying to minimize the duration of the trajectory,
	% so the cost is proportional to the duration.
	% traj_create_cost() takes in the (symbolic) cost function. It also allows a name to be specified
	% for the cost.
	cost = traj_create_cost(duration, 'Time');
	scenario = traj_add_cost(scenario, cost);

	% Run the optimization
	% This returns the scenario and a success flag (also stored in the scenario)
	% that is true if it was successful and false otherwise.
	[scenario,success] = traj_optimize(scenario);

	% Check if it failed -- if so, run the constraint feasibility analyzer
	if ~success
		% It failed -- run the analyzer
		traj_analyze_feasibility(scenario);
	end
end
