% This function clones the given phase. You pass it the original phase and the name for the new phase

function phase2 = traj_clone_phase(phase1, name)
	disp(['Cloning phase ''' phase1.names.phase ''' into phase ''' name '''']);

	% For now, this is just a simply copy operation...
	phase2             = phase1;
	phase2.names.phase = name;
end
