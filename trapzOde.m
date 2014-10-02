% This function implements the Trapezoidal method for solving ODES.
% It adds a constraint to the given problem to enforce the derivative/integral relationship
%
% Parameters:
%     prob The OptTool object to add the constraint to
%     y    The variable whose trajectory is defined by the ODE
%     dy   The derivative of y
%     T    Total duration of the trajectory

function trapzOde(prob, y, dy, T)
	% Calculate the time between nodes
	dt = T / (numel(y)-1);

	% The trapezoidal method defines the differences between y's nodes
	% by using averaged values from dy.
	delta_y = y(2:end) - y(1:end-1);

	pred_delta_y = dt * (dy(1:end-1) + dy(2:end))/2;

	% Last, add in the constraint, setting actual and predicted deltas equal to each other
	prob.addCon(delta_y, '==', pred_delta_y)
end
