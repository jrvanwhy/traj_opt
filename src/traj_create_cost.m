% This function creates a cost structure. Parameters:
%     name  The name of this cost
%     value The value of the cost

% Structure of the cost struct:
% name        The cost's name
% fcn         The cost function
% struct_type The type of the struct. Always 'cost'
% version     The version of this struct.

% Whether fcn is an anonymous function or a symbolic expression depends on the context
% in which the cost resides. As returned by the user's functions, it's a symbolic
% expression -- but afterwards, it is converted to a function handle of parameters that depend
% on its context (defined by the larger structure it is encapsulated within).

function cost = traj_create_cost(name, value)
	% Generate struct type and version
	cost.version     = traj_version();
	cost.struct_type = 'cost';

	% Copy over the name and value (simplifying the value)
	cost.name = name;
	cost.fcn  = simplify(value);
end
