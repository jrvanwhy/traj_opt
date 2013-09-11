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
% Also, the number of "nodes" for the dynamics representation (equivalent to timesteps for an ODE solver)
% is a parameter of this function.
function scenario = car_time(n_nodes)
	% Create the scenario structure itself. This doesn't require any parameters.
	scenario = traj_create_scenario();

	% This creates a state structure. We specify names for each element in the state.
	% The state will be a column vector, even if the names array isn't. The
	% number of elements in the names array determines the size of the state.
	% We're using a cell array here because these are strings
	% If a string is passed in, it the state should be of size 1
	state = traj_create_state({ 'position', 'velocity' });

	% This is a similar to traj_create_state(), but for input instead.
	input = traj_create_input('force');

	% This is our dynamics function. It takes in the state and input and returns the derivative
	% of the state. It gets called by the optimization framework (during the call to traj_create_dynamics())
	function dx = dynamics(x, u)
		% We'll use a linear representation of the system here
		A = [ 0 1
		      0 0 ];

		B = [ 0
		      1 ];

		% Specify the derivative of the state
		dx = A * x + B * u;
	end

	% Create the dynamics structure
	dynamics = traj_create_dynamics(@dynamics, state, input);

	% Here, we create a dynamic phase. We must specify the dynamics
	% in use and may specify the number of nodes for this phase.
	% Additionally, we may create a "cost" function via traj_create_cost()
	% or constraints and add them to the resulting dynamic phase.
	phase = traj_create_phase(dynamics, n_nodes);

	% Now we need to constrain the starting and ending states to be equal to the origin
	% and one unit in the positive direction, respectively. To do this, we first need to be able
	% to access the states. We can get a vector of all the states using traj_get_states(phase)
	states = traj_get_states(phase);

	% Create the constraint for the starting state. Set it to the origin. traj_create_constraint()
	% takes three parameters. The 2nd parameter sets the constraint type ('='/'==', '<=', or '>='), and the other
	% two are the "sides" of the equation or inequality.
	initial_state_constraint = traj_create_constraint(states(:,1), '=', [ 0
	                                                                      0 ]);

	% Similarly, create an inequality constraint for the ending state
	ending_state_constraint = traj_create_constraint(states(:,end), '.', [ 1
	                                                                       0 ]);

	% Add the constraints to the phase
	% traj_add_constraint() adds the constraint to the given phase or scenario
	phase = traj_add_constraint(phase, initial_state_constraint);
	phase = traj_add_constraint(phase, ending_state_constraint);

	% Add the phase to the scenario
	scenario = traj_add_phase(scenario, phase);

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
