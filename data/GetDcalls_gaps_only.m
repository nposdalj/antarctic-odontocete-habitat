%% Recover missing detections using OLD gaps
% Workflow:
% 1. Concatenate all OLD detections
% 2. Extract NEW file start times from filenames
% 3. Identify gaps >30 mins before each NEW file start
% 4. Keep ONLY NEW detections inside those gaps
% 5. Save filtered NEW .mat files with SAME filenames

clear; clc;

%% --------------------------------------------------------
%% USER SETTINGS
%% --------------------------------------------------------

site = 'EIE';

% OLD detector MAT files
oldDir = ['\\snowman.ucsd.edu\Ally_Working_Disk\Analysis\Bm\Bm D call detector output\Antarc\',site,'\data'];

% NEW detector MAT files
newDir = ['F:\Antarctic_Mysticetes\D calls\',site,'\filtered'];

% Output folder for filtered MAT files
outDir = ['F:\Antarctic_Mysticetes\D calls\',site,'_gaps_only'];

if ~exist(outDir,'dir')
    mkdir(outDir)
end

% Final recording end time
recordingEnd = datetime( ...
    '2016-02-15 00:00:00', ...
    'TimeZone','UTC');

% Gap threshold
gapThresh = minutes(30);

%% --------------------------------------------------------
%% SITE RECORDING PERIODS
%% --------------------------------------------------------

switch site

    case 'EIE'

        recordingStart = datetime( ...
            '2016-02-03 23:59:59', ...
            'TimeZone','UTC');

        recordingEnd = datetime( ...
            '2016-12-02 05:15:05', ...
            'TimeZone','UTC');

    case 'SSI'

        recordingStart = datetime( ...
            '2015-02-10 00:00:00', ...
            'TimeZone','UTC');

        recordingEnd = datetime( ...
            '2016-01-29 13:15:25', ...
            'TimeZone','UTC');

    case 'EI'

        recordingStart = datetime( ...
            '2014-03-05 00:00:00', ...
            'TimeZone','UTC');

        recordingEnd = datetime( ...
            '2014-07-14 12:03:50', ...
            'TimeZone','UTC');

    otherwise

        error('Unknown site');

end

%% --------------------------------------------------------
%% LOAD + CONCATENATE OLD DETECTIONS
%% --------------------------------------------------------

fprintf('\n====================================\n');
fprintf('LOADING OLD DETECTIONS\n');
fprintf('====================================\n');

oldFiles = dir(fullfile(oldDir,'*.mat'));

allOldTimes = datetime.empty(0,1);
allOldTimes.TimeZone = 'UTC';

for k = 1:numel(oldFiles)

    try

        fprintf('Loading %s\n', oldFiles(k).name);

        fpath = fullfile(oldFiles(k).folder, ...
            oldFiles(k).name);

        load(fpath,'-mat');

        calls = hyd.detection.calls;

        [tStartRaw,~] = getStartEnd(calls);

        tStart = toDatetime(tStartRaw);

        allOldTimes = [allOldTimes; tStart];

    catch ME

        warning('Problem loading OLD file:\n%s', ...
            ME.message);

    end
end

%% Sort OLD detections

allOldTimes = sort(allOldTimes);
%% Restrict OLD detections to deployment period

idxValid = allOldTimes >= recordingStart & ...
           allOldTimes <= recordingEnd;

allOldTimes = allOldTimes(idxValid);

fprintf('\nTotal OLD detections: %d\n', ...
    numel(allOldTimes));

%% --------------------------------------------------------
%% GET NEW FILE START TIMES
%% --------------------------------------------------------

fprintf('\n====================================\n');
fprintf('READING NEW FILE START TIMES\n');
fprintf('====================================\n');

newFiles = dir(fullfile(newDir,'*.mat'));

newStarts = datetime.empty(0,1);
newStarts.TimeZone = 'UTC';

for k = 1:numel(newFiles)

    fname = newFiles(k).name;

    token = regexp(fname, ...
        '_(\d{6})_(\d{6})\.d\d+\.mat$', ...
        'tokens');

    if isempty(token)

        warning('Could not parse filename:\n%s', fname);

        continue

    end

    token = token{1};

    startStr = [token{1} token{2}];

    dt = datetime( ...
        startStr, ...
        'InputFormat','yyMMddHHmmss', ...
        'TimeZone','UTC');

    newStarts(end+1,1) = dt;

end

%% Sort NEW starts

[newStarts,idx] = sort(newStarts);

newFiles = newFiles(idx);

%% --------------------------------------------------------
%% FIND GAPS
%% --------------------------------------------------------

fprintf('\n====================================\n');
fprintf('IDENTIFYING GAPS\n');
fprintf('====================================\n');

gapStarts = datetime.empty(0,1);
gapStarts.TimeZone = 'UTC';

gapEnds = datetime.empty(0,1);
gapEnds.TimeZone = 'UTC';

%% ----- Gaps before each NEW file start -----

for k = 2:numel(newStarts)

    thisStart = newStarts(k);

    % Find OLD detections before this start
    idxOld = allOldTimes < thisStart;

    if ~any(idxOld)
        continue
    end

    prevOld = max(allOldTimes(idxOld));

    gapDur = thisStart - prevOld;

    if gapDur > gapThresh

        fprintf('\nGAP FOUND\n');

        fprintf('Gap start: %s\n', ...
            datestr(prevOld));

        fprintf('Gap end:   %s\n', ...
            datestr(thisStart));

        fprintf('Duration: %.2f hrs\n', ...
            hours(gapDur));

        gapStarts(end+1,1) = prevOld;
        gapEnds(end+1,1)   = thisStart;

    end
end

%% ----- Final recording gap -----

idxLast = allOldTimes < recordingEnd;

if any(idxLast)

    lastOld = max(allOldTimes(idxLast));

    finalGap = recordingEnd - lastOld;

    if finalGap > gapThresh

        fprintf('\nFINAL GAP FOUND\n');

        fprintf('Gap start: %s\n', ...
            datestr(lastOld));

        fprintf('Gap end:   %s\n', ...
            datestr(recordingEnd));

        fprintf('Duration: %.2f hrs\n', ...
            hours(finalGap));

        gapStarts(end+1,1) = lastOld;
        gapEnds(end+1,1)   = recordingEnd;

    end
end

%% --------------------------------------------------------
%% SAVE GAP TABLE
%% --------------------------------------------------------

GapTable = table( ...
    gapStarts, ...
    gapEnds, ...
    gapEnds-gapStarts, ...
    'VariableNames', ...
    {'GapStart','GapEnd','Duration'});

gapCSV = fullfile(outDir, ...
    [site '_GapSummary.csv']);

writetable(GapTable,gapCSV);

fprintf('\nSaved gap summary:\n%s\n', gapCSV);

%% --------------------------------------------------------
%% FILTER NEW FILES
%% --------------------------------------------------------

fprintf('\n====================================\n');
fprintf('FILTERING NEW DETECTIONS\n');
fprintf('====================================\n');

for k = 1:numel(newFiles)

    try

        fprintf('\nProcessing %s\n', ...
            newFiles(k).name);

        fpath = fullfile(newFiles(k).folder, ...
            newFiles(k).name);

        S = load(fpath,'-mat');

        calls = S.hyd.detection.calls;

        [tStartRaw,~] = getStartEnd(calls);

        detTimes = toDatetime(tStartRaw);

        %% ----- KEEP ONLY DETECTIONS INSIDE GAPS -----

% First remove detections outside deployment period

keepIdx = detTimes >= recordingStart & ...
          detTimes <= recordingEnd;

% Then apply gap filtering

gapKeep = false(size(detTimes));
        for g = 1:numel(gapStarts)

    gapKeep = gapKeep | ...
    (detTimes >= gapStarts(g) & ...
                 detTimes <= gapEnds(g));

        end

        keepIdx = keepIdx & gapKeep;

        fprintf('Keeping %d detections\n', ...
            sum(keepIdx));

        %% ----- FILTER CALLS -----

        if istable(calls)

            callsFiltered = calls(keepIdx,:);

        elseif isstruct(calls)

            if numel(calls) > 1

                callsFiltered = calls(keepIdx);

            else

                fields = fieldnames(calls);

                for f = 1:numel(fields)

                    thisField = calls.(fields{f});

                    if numel(thisField) == numel(keepIdx)

                        callsFiltered.(fields{f}) = ...
                            thisField(keepIdx);

                    else

                        callsFiltered.(fields{f}) = ...
                            thisField;

                    end
                end
            end
        end

        %% ----- SAVE SAME FILENAME -----

        S.hyd.detection.calls = callsFiltered;

        outFile = fullfile(outDir, ...
            newFiles(k).name);

        save(outFile,'-struct','S');

        fprintf('Saved:\n%s\n', outFile);

    catch ME

        warning('Problem processing NEW file:\n%s', ...
            ME.message);

    end
end

fprintf('\n====================================\n');
fprintf('FINISHED\n');
fprintf('====================================\n');

%% ========================================================
%% HELPER FUNCTIONS
%% ========================================================

function [tStartRaw,tEndRaw] = getStartEnd(calls)

    if istable(calls)

        tStartRaw = calls.julian_start_time;
        tEndRaw   = calls.julian_end_time;

        return
    end

    if isstruct(calls)

        if numel(calls) > 1

            tStartRaw = [calls.julian_start_time]';
            tEndRaw   = [calls.julian_end_time]';

        else

            tStartRaw = calls.julian_start_time;
            tEndRaw   = calls.julian_end_time;

        end

        return
    end

    error('Unsupported calls type');

end

%% --------------------------------------------------------

function dt = toDatetime(x)

    if isempty(x)

        dt = datetime.empty(0,1);
        dt.TimeZone = 'UTC';

        return
    end

    if isdatetime(x)

        dt = x;

        if isempty(dt.TimeZone)
            dt.TimeZone = 'UTC';
        end

        return
    end

    if isnumeric(x)

        x = x(:);

        xmax = max(x);

        if xmax > 1e9

            dt = datetime( ...
                x, ...
                'ConvertFrom','posixtime', ...
                'TimeZone','UTC');

        elseif xmax > 5e5

            dt = datetime( ...
                x, ...
                'ConvertFrom','datenum', ...
                'TimeZone','UTC');

        else

            error('Unknown numeric time format');

        end

        return
    end

    error('Unsupported time type');

end