% ==============================
% USER INPUT
% ==============================
data_dir = 'F:\Antarctic_Mysticetes\D calls\filtered';  % <-- your folder

files = dir(fullfile(data_dir, '*.mat'));

all_t = [];  % store all detection times

% ==============================
% LOAD ALL FILES
% ==============================
for i = 1:length(files)
    
    fprintf('Loading %d/%d: %s\n', i, length(files), files(i).name);
    
    S = load(fullfile(data_dir, files(i).name));
    
    if ~isfield(S, 'hyd') || isempty(S.hyd)
        continue
    end
    
    if ~isfield(S.hyd, 'detection') || isempty(S.hyd.detection)
        continue
    end
    
    calls = S.hyd.detection.calls;
    
    if isempty(calls)
        continue
    end
    
    % Extract times
    t = [calls.julian_start_time];
    
    all_t = [all_t, t]; %#ok<AGROW>
end

% ==============================
% SORT + CONVERT
% ==============================
all_t = sort(all_t);
t_dt = datetime(all_t, 'ConvertFrom', 'datenum');

% ==============================
% PLOT
% ==============================
figure;
plot(t_dt, ones(size(t_dt)), '.k')
ylim([0.9 1.1])
yticks([])

xlabel('Date')
title('All Detection Timeline')
grid on

dt = diff(t_dt);

gap_idx = find(dt > days(1));  % gaps > 1 day

fprintf('\nFound %d large gaps:\n', length(gap_idx));

for i = 1:length(gap_idx)
    fprintf('Gap from %s to %s (%.2f days)\n', ...
        datestr(t_dt(gap_idx(i))), ...
        datestr(t_dt(gap_idx(i)+1)), ...
        days(dt(gap_idx(i))));
end
