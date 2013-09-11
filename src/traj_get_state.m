% This function retrieves the state (or a subset of the state) at the given time(s).
% If substates is not passed, it returns the full state

function state = traj_get_state(scenario, t, substates)
	% Set substates to the whole state, if it was not passed in.
	if nargin < 3
		substates = 1:size(scenario.states.def_state, 1);
	end

	% Each optimization method requires a different function to retrieve the state at a given time.
	% Switch based on the method, error if the method isn't recognized.
	switch scenario.method
		case 'direct collocation'; state = traj_get_state_dircol(scenario, t, substates);
		otherwise; error(['scenario.method of ', scenario.method, ' not recognized!']);
	end
end

% This gets the state at a given time for direct collocation, using the cubic polynomial approximation to the state
function state = traj_get_state_dircol(scenario, t, substates)
	% Identify which shoot this is, including a correction for the endpoints.
	dt       = scenario.num_duration / scenario.shoots;  % Time per shoot
	shoot    = floor(t / dt) + 1;                        % Shoot number. This may exceed the total number of shoots (if t = duration)
	shoot    = min(shoot, scenario.shoots);              % Now it should always be valid.

	% Time relative to the start of the shoot
	rel_t    = t - dt * (shoot - 1);

	% State at the beginning and end of the shoot, as well as the derivative of the state
	% at the beginning and end of the shoot. These only include the relative substates.
	x0  = scenario.states.num_states( substates, 2 * shoot - 1);
	x1  = scenario.states.num_states( substates, 2 * shoot + 1);
	dx0 = scenario.states.num_dstates(substates, 2 * shoot - 1);
	dx1 = scenario.states.num_dstates(substates, 2 * shoot + 1);

	% Do the polynomial interpolation
	% This approximation was obtained by solving the following equations for a, b, c, and d:
	% p(t) = a*t^3 + b*t^2 + c*t + d
	% ṗ(t) = 3*a*t^2 + 2*b*t + c
	% p(0)  = x0
	% p(dt) = x1
	% ṗ(0)  = dx0
	% ṗ(dt) = dx1
	b = (3 * (x1 - x0) / dt - dx1 - 2*dx0) / dt;
	c = dx0;
	d = x0;
	a = (((x1 - x0)/dt - dx0)/dt - b)/dt;
	state = ((a*rel_t + b)*rel_t + c)*rel_t + d; % a*t^3 + b*t^2 + c*t + d == ((a*t + b)*t + c)*t + d
end
