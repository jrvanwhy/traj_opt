% This function does most of the setup and processing for a given scenario.
% It returns a scenario ready for optimization (or generation
% of a fast optimization function, if desired).

% Scenario structure setup:
% name        This scenario's name. Blank if not specified by the user.
% version     The version of the optimizer that generated this scenario
% struct_type The type of this structure -- always 'scenario'

function scenario = traj_setup_scenario(varargin)
	disp('Initializing scenario structure')

	% Basic structure initialization
	scenario.version     = traj_version();
	scenario.struct_type = 'scenario';

	% Initialize the name to a blank string -- it will be changed later
	% if the user specified a name
	scenario.name = '';

	% Process each input in order
	for iter = 1:numel(varargin)
		% Switch based on the type of this input
		if ischar(varargin{iter})
			% This is a name for the scenario
			disp(['	Processing name ''' varargin{iter} ''''])

			% Check if a name's already been supplied -- if so, error
			if ~isempty(scenario.name)
				error(['Multiple names specified. The first: ''' scenario.name ''', the second: ''' varargin{iter} ''''])
			end

			% Add the name to the scenario
			scenario.name = varargin{iter};
		elseif isstruct(varargin{iter}) && isfield(varargin{iter}, 'struct_type') && varargin{iter}.struct_type == 'dynamic phase'
			% This is a new phase.
			disp(['	Processing phase ''' varargin{iter}.name ''''])

			% Call the phase processing function
			scenario = add_phase(scenario, varargin{iter});
		else
			% Spit out an error -- this input is not recognized
			error(['Input ' num2str(iter) ' not recognized'])
		end
	end
end

% This function adds a phase to the given scenario
function scenario = add_phase(scenario, phase)
	error('TODO: this')
end
