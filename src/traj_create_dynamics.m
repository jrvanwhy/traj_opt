% This function initializes a new dynamics structure. Its arguments are a function handle to the dynamics function,
% the state associated with this set of dynamics, and the input associated with this set of dynamics.

% The dynamics structure is documented in documentation.pdf in the 'doc' directory.

function dynamics = traj_create_dynamics(dynamicsfcn, state, input)
	% Check if the proper number of arguments was specified
	if nargin < 3
		error(['3 parameters are needed, only ', num2str(nargin), ' passed.'])
	end

	% Verify dynamicsfcn is a function handle
	if ~isa(dynamicsfcn, 'function_handle')
		error('dynamicsfcn argument must be a function handle.')
	end

	% Verify state is a state structure
	if ~isstruct(state) || ~isfield(state, 'struct_type') || ~strcmp(state.struct_type, 'state')
		error('state must be a state structure, as created by traj_create_state()')
	end

	% Verify input is an input structure
	if ~isstruct(input) || ~isfield(input, 'struct_type') || ~strcmp(input.struct_type, 'input')
		error('input must be a input structure, as created by traj_create_input()')
	end

	% Fill in the version and struct type
	dynamics.version     = traj_version();
	dynamics.struct_type = 'dynamics';

	% Simply assign the state and input fields
	dynamics.state = state;
	dynamics.input = input;

	% Computing the value of the dynamics_fcn field is a multistep process. First, we obtain the
	% symbolic form of the state and input. Next we call dynamicsfcn with the given state and input, then simplify it.
	% Last, we use matlabFunction to make it into an anonymous function (using matlabFunction to create a function then
	% calling that function is still faster than using subs() on the expression, in my testing).

	% Create symbolic representations of the state and input
	sym_state = sym('ex_x', size(state.names));
	sym_input = sym('ex_u', size(input.names));

	% Call the given dynamics function with the example state and input
	dynamics_expr = dynamicsfcn(sym_state, sym_input);

	% Simplify the dynamics equations
	dynamics_expr = simplify(dynamics_expr);

	% Convert dynamics_expr into a function of state and input
	dynamics.dynamics_fcn = matlabFunction(dynamics_expr, 'vars', {sym_state, sym_input});

	% Clean up the symbolic variables
	syms sym_state sym_input clear
end
