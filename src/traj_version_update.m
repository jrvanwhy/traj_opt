% This function checks the version of the structure, and updates it if necessary.

function traj_struct = traj_version_update(traj_struct)
	% I haven't implemented this yet, since it's still on version 1. Check if the structure version
	% doesn't match the current version -- if so, error.
	if traj_struct.version ~= traj_version()
		error('TODO: this')
	end

	% Make MATLAB happy
	traj_struct = traj_struct;
end
