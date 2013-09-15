% This function initializes a new input structure. Its arguments are the names of the
% sub-inputs of this input.

% The input structure is documented in documentation.pdf in the 'doc' directory.

% We use varargin so the user may specify any natural number of parameters (at their choosing).
% This is MATLAB's form of variadic functions.
function input = traj_create_input(varargin)
	% Error checking on input -- there must be at least one input, and all parameters must be strings.
	if size(varargin, 2) == 0
		error('At least one input variable must be created. No names were given!');
		return;
	end

	% Check if all parameters are strings
	if ~iscellstr(varargin)
		error('All parameters must be strings.');
		return;
	end

	% Fill in the version and struct type
	input.version     = traj_version();
	input.struct_type = 'input';

	% Fill in the names cell array. Since varargin is 1 by N but names is N by 1, we must transpose
	% the varargin cell array.
	input.names = varargin.';
end
