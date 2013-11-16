% This function creates a constraint structure. Parameters:
%     name  The name of this constraint
%     lside The left hand side of the constraint equation or inequality
%     type  The constraint's 'type' -- '<=', '=', '==', or '>='
%     rside The right hand side of the constraint equation or inequality

% Structure of the constraint struct:
% name        The constraint's name
% fcn         The constraint function (g(x) <= 0 or h(x) = 0)
% struct_type The type of the struct. Always 'constraint'
% con_type    The type of the constraint ('<=' or '=')
% version     The version of this struct.

% Whether fcn is an anonymous function or a symbolic expression depends on the context
% in which the constraint resides. As returned by the user's functions, it's a symbolic
% expression -- but afterwards, it is converted to a function handle of parameters that depend
% on its context (defined by the larger structure it is encapsulated within).

function constraint = traj_create_constraint(name, lside, type, rside)
	% Error checking
	if size(lside) ~= size(rside)
		error('lside and rside of constraint must have equal sizes')
	end

	% Generate name and version
	constraint.version     = traj_version();
	constraint.struct_type = 'constraint';

	% Copy over the name
	constraint.name = name;

	% Convert the two sides of the constraint into a single function through
	% additive inversion and summation (subtract rside from lside ;p)
	constraint.fcn = lside - rside;

	% Catch the case of a constraint of the form '>=' -- in this case, we need to
	% negate the constraint's function and set the type back to '<=' to match
	% the typical optimization problem constraint setup.
	if type == '>='
		constraint.fcn = -constraint.fcn;
		type = '<=';
	end

	% Alias '==' to '='
	if type == '=='
		type = '=';
	end

	% Simplify the constraint function
	constraint.fcn = simplify(constraint.fcn);

	% Copy over the constraint's type
	constraint.con_type = type;
end
