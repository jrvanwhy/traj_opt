% This function initializes a new state structure. Its arguments are the names of the
% sub-states of this state.

% The state structure is documented in documentation.pdf in the 'doc' directory.

% We use varargin so the user may specify any natural number of parameters (at their choosing).
% This is MATLAB's form of variadic functions.
function state = traj_create_state(varargin)
	% Error checking on input -- there must be at least one input, and they must all be strings.
	if size(varargin, 2) == 0
		error('At least one state variable must be created. No names were given!');
		return;
	end

	% Check if all inputs are strings
	if ~iscellstr(varargin)
		error('All inputs must be strings.');
		return;
	end

	% Fill in the version and struct type
	state.version     = traj_version();
	state.struct_type = 'state';

	% Fill in the names cell array. Since varargin is 1 by N but names is N by 1, we must transpose
	% the varargin cell array.
	state.names = varargin.';
end
