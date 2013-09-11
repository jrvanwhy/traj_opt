% This function initializes the field of a structure to value if it has not been declared (or empty if value was not given)
% If it has been declared already, it does nothing

function struct = declare_field(struct, fieldname, value)
	if isfield(struct, fieldname)
		% We have to set struct here or it'll give us an error.
		struct = struct;
	else
		% Dynamically add the new field; initialize as value (or empty, if value is not given).
		if nargin < 3
			value = [];
		end
		struct = setfield(struct, fieldname, value);
	end
end
