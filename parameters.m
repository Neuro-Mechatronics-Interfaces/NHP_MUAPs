function varargout = parameters(varargin)
%PARAMETERS Return parameters struct, which sets default values for things like epoch durations etc.
%
% Example 1
%   pars = parameters(); % Returns full "default parameters" struct.
%
% Example 2
%   f = parameters('raw_data_folder'); % Returns raw data folder value
%
% I tried to pull the key parameters, which are really path information
% from your local mapping to where the raw data is and where any
% auto-exported figures/data objects should go. Both
% `generated_data_folder` and `raw_data_folder` should be the folder that
% contains "animal name" folders (e.g. the folder with `Forrest` in it). 
%
% See also: Contents, load_tmsi_raw

pars = struct;
pars.data_share_root = 'R:/NMLShare';
pars.raw_data_folder = fullfile(pars.data_share_root, 'raw_data/primate');
pars.generated_data_folder = fullfile(pars.data_share_root, 'generated_data/human/TC');
pars.task_bits_file = "wrist_task_bits.json";pars.version = 1.0; % "Version" of repo code -- check that this matches the TAG

N = numel(varargin);
if nargout == 1
    if rem(N, 2) == 1
        varargout = {pars.(varargin{end})};
        return;
    else
        f = fieldnames(pars);
        for iV = 1:2:N
            idx = strcmpi(f, varargin{iV});
            if sum(idx) == 1
               pars.(f{idx}) = varargin{iV+1}; 
            end
        end
        varargout = {pars};
        return;
    end
else
    f = fieldnames(pars);
    varargout = cell(1, nargout);
    for iV = 1:numel(varargout)
        idx = strcmpi(f, varargin{iV});
        if sum(idx) == 1
            varargout{iV} = pars.(f{idx}); 
        else
            error('Could not find parameter: %s', varargin{iV}); 
        end
    end
end

end
