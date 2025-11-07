%% Batch extract Blue Whale D-call times from .mat -> CSV
% Configure paths
site = 'SSI';
inDir  = ['\\snowman.ucsd.edu\Ally_Working_Disk\Analysis\Bm\Bm D call detector output\Antarc\',site,'\Data'];
outDir = 'L:\Shared drives\Antarctic Marine Mammals\Marine Mammal Data\Mysticetes';
outCSV = fullfile(outDir, [site,'_BlueWhale_D_calls.csv']);

% Scan files
files = dir(fullfile(inDir, '*.mat'));
allStart = datetime.empty(0,1);  allStart.TimeZone = 'UTC';
allEnd   = datetime.empty(0,1);  allEnd.TimeZone   = 'UTC';
allSrc   = strings(0,1);

for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    load(fpath, '-mat');
    calls = hyd.detection.calls;

    [tStartRaw, tEndRaw] = getStartEnd(calls); % handles table/struct/struct-array
    tStart = toDatetime(tStartRaw);
    tEnd   = toDatetime(tEndRaw);

    % Ensure column vectors and matching lengths
    n = min(numel(tStart), numel(tEnd));
    tStart = tStart(1:n);
    tEnd   = tEnd(1:n);

    allStart = [allStart; tStart]; %#ok<AGROW>
    allEnd   = [allEnd;   tEnd];   %#ok<AGROW>
    allSrc   = [allSrc;   repmat(string(files(k).name), n, 1)]; %#ok<AGROW>
end

% Assemble, sort, and export
T = table(allStart, allEnd, 'VariableNames', {'StartTime','EndTime'});
T = sortrows(T, 'StartTime');
writetable(T, outCSV);
fprintf('Wrote %d calls to %s\n', height(T), outCSV);

%% ---------- Helpers Functions ----------
function [tStartRaw, tEndRaw] = getStartEnd(calls)
% Return raw julian_start_time / julian_end_time from a table or struct/struct-array.
    if istable(calls)
        assert(any(strcmpi(calls.Properties.VariableNames,'julian_start_time')) && ...
               any(strcmpi(calls.Properties.VariableNames,'julian_end_time')), ...
               'Calls table lacks julian_start_time or julian_end_time.');
        tStartRaw = calls.(calls.Properties.VariableNames{strcmpi(calls.Properties.VariableNames,'julian_start_time')});
        tEndRaw   = calls.(calls.Properties.VariableNames{strcmpi(calls.Properties.VariableNames,'julian_end_time')});
        return
    end

    if isstruct(calls)
        assert(isfield(calls,'julian_start_time') && isfield(calls,'julian_end_time'), ...
               'Calls struct lacks julian_start_time or julian_end_time.');
        if numel(calls) > 1
            % struct array with scalar fields
            tStartRaw = [calls.julian_start_time]';
            tEndRaw   = [calls.julian_end_time]';
        else
            % scalar struct with vector fields
            tStartRaw = calls.julian_start_time;
            tEndRaw   = calls.julian_end_time;
        end
        return
    end

    error('Unsupported calls type: %s', class(calls));
end

function dt = toDatetime(x)
% Convert datetime / datenum / POSIX seconds to datetime (UTC)
    if isempty(x); dt = []; return; end
    if isdatetime(x)
        dt = x;
        if isempty(dt.TimeZone); dt.TimeZone = 'UTC'; end
        return
    end
    if isnumeric(x)
        x = x(:);
        xmax = max(x);
        if xmax > 1e9
            dt = datetime(x, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');   % POSIX seconds
        elseif xmax > 5e5
            dt = datetime(x, 'ConvertFrom', 'datenum',   'TimeZone', 'UTC');   % MATLAB datenum
        else
            error('Values look like samples/seconds offset. Convert using recording start and/or sample rate first.');
        end
        return
    end
    error('Unsupported time type: %s', class(x));
end