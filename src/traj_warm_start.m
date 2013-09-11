% This function sets up the scenario to be optimized again (perhaps with different additional parameters)

function scenario = traj_warm_start(scenario)
	scenario.params.init_params = scenario.minimum(:,end);
end
