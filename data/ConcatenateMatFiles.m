% ==============================
% USER INPUT
% ==============================
data_dir = 'F:\Antarctic_Mysticetes\D calls\EIE\filtered';  % your folder
out_name = 'All_Detections_Concatenated';

files = dir(fullfile(data_dir, '*.mat'));

all_start = [];
all_end   = [];

% ==============================
% LOAD + CONCATENATE
% ==============================
for i = 1:length(files)
    
    fprintf('Loading %d/%d: %s\n', i, length(files), files(i).name);
    
    S = load(fullfile(data_dir, files(i).name));
    
    if ~isfield(S, 'hyd') || isempty(S.hyd)
        continue
    end
    
    calls = S.hyd.detection.calls;
    
    if isempty(calls)
        continue
    end
    
    % Extract times
    start_t = [calls.julian_start_time];
    end_t   = [calls.julian_end_time];
    
    all_start = [all_start, start_t]; %#ok<AGROW>
    all_end   = [all_end, end_t]; %#ok<AGROW>
end

% ==============================
% SORT BY TIME
% ==============================
[all_start, sort_idx] = sort(all_start);
all_end = all_end(sort_idx);

% ==============================
% CONVERT TO DATETIME
% ==============================
StartTime = datetime(all_start, 'ConvertFrom', 'datenum');
EndTime   = datetime(all_end,   'ConvertFrom', 'datenum');

% ==============================
% BUILD TABLE (MATCH CSV STYLE)
% ==============================
T = table(StartTime(:), EndTime(:));

% ==============================
% SAVE OUTPUTS
% ==============================
save(fullfile(data_dir, [out_name '.mat']), 'T')
writetable(T, fullfile(data_dir, [out_name '.csv']))

fprintf('\nSaved concatenated outputs:\n');
fprintf('%s.mat\n', out_name);
fprintf('%s.csv\n', out_name);