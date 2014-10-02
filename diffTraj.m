% Copyright (C) 2014 Johnathan Van Why
% See LICENSE.txt for details

% This function differentiates a trajectory as a function of time.
% It optionally creates second derivatives as well.
%
% Parameters:
%     y The trajectory to be differentiated
%     T The overall trajectory duration
%
% Returns:
%     dy  The derivative of y as a function of time
%     ddy The second derivative of y as a function of time

function [dy,ddy] = diffTraj(y, T)
	% Calculate the delta time between the nodes
	dt = T / (numel(y)-1);

	% Do central differencing (except at the endpoints, where forward and backwards finite
	% differences are required)
	dy = [y(2)-y(1); (y(3:end) - y(1:end-2))/2; y(end)-y(end-1)]/dt;

	% Exit if one derivative was all that was required
	if nargout < 2
		return
	end

	% Do second-order differentiation as well. This time, duplicate the
	% second and second-to-last values to handle the endpoints.
	ddy = (y(1:end-2) - 2*y(2:end-1) + y(3:end))/dt^2;
	ddy = [ddy(1); ddy; ddy(end)];
end
