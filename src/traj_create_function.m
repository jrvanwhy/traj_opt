% This function creates a function structure. Parameters:
%     name  The name of this function
%     value The value of the function

% Structure of the function struct:
% name        The function's name
% fcn         The function function
% struct_type The type of the struct. Always 'function'
% version     The version of this struct.

% Whether fcn is an anonymous function or a symbolic expression depends on the context
% in which the function resides. As returned by the user's functions, it's a symbolic
% expression -- but afterwards, it is converted to a function handle of parameters that depend
% on its context (defined by the larger structure it is encapsulated within).

function func = traj_create_function(name, value)
	% Generate struct type and version
	func.version     = traj_version();
	func.struct_type = 'function';

	% Copy over the name and value (simplifying the value)
	func.name = name;
	func.fcn  = simplify(value);
end
