% This function creates the scenario structure. If you pass it phases or constraints, it will initialize the scenario
% with those phases and constraints.

% The scenario, phase, and constraint structures are documented in documentation.pdf in the 'doc' directory.

% We use varargin so the user may specify any natural number of parameters (at their choosing).
% This is MATLAB's form of variadic functions.
function scenario = traj_create_scenario(varargin)
	% Fill in the version field of the struct.
	scenario.version = traj_version();

	% Fill in the struct_type field, as spec'd in the documentation for the scenario structure.
	scenario.struct_type = 'scenario';

	% At this point, we should go through every argument passed in, update it (via the versioning
	% system I haven't written yet), then add it to scenario using the relevant function.
	% Since this hasn't been written yet, we'll choose to check for arguments, and if any were passed,
	% spit out an error and die.
	if size(varargin, 2) > 0
		error('traj_create_scenario(): Adding constraints/phases to the scenario during creation not yet possible!');
	end
end
