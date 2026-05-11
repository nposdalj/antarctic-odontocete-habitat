% =========================================================
% USER INPUT
% =========================================================
xwav_dir = 'U:\Antarctica\Antarc_SSI_01\Antarc_SSI_01_disks01-08_df100';   % <-- change this
n_files = 5;                      % number of files to check

% =========================================================
% GET FILE LIST
% =========================================================
files = dir(fullfile(xwav_dir, '*.x.wav'));
files = sort({files.name});

n_files = min(n_files, length(files));

fprintf('Checking first %d XWAV files...\n\n', n_files);

% =========================================================
% STORAGE
% =========================================================
start_times = datetime.empty(n_files,0);
durations_sec = nan(n_files,1);

% =========================================================
% LOOP THROUGH FILES
% =========================================================
for i = 1:n_files
    fname = fullfile(xwav_dir, files{i});
    
    % Read audio info
    info = audioinfo(fname);
    
    % Extract start time from filename OR header
    % ---- OPTION 1: parse from filename (common format) ----
    % Example: 20150715_123000.x.wav
    tokens = regexp(files{i}, '(\d{6})_(\d{6})', 'tokens');    
    if ~isempty(tokens)
        date_str = tokens{1}{1};
        time_str = tokens{1}{2};
        start_times(i) = datetime([date_str time_str], ...
            'InputFormat','yyyyMMddHHmmss');
    else
        % ---- OPTION 2: fallback (if Triton header needed) ----
        warning('Could not parse time from filename: %s', files{i});
        start_times(i) = NaT;
    end
    
    % Duration from samples
    durations_sec(i) = info.Duration;
    
    fprintf('File %d: %s\n', i, files{i});
    fprintf('  Start time: %s\n', datestr(start_times(i)));
    fprintf('  Duration: %.2f sec\n\n', durations_sec(i));
end

% =========================================================
% COMPUTE DELTAS
% =========================================================
fprintf('\n==============================\n');
fprintf('TIME STEP CHECK\n');
fprintf('==============================\n');

for i = 1:n_files-1
    dt = seconds(start_times(i+1) - start_times(i));
    
    fprintf('Between file %d -> %d:\n', i, i+1);
    fprintf('  Delta t: %.2f sec\n', dt);
    fprintf('  Expected: %.2f sec\n', durations_sec(i));
    fprintf('  Difference: %.2f sec\n\n', dt - durations_sec(i));
end