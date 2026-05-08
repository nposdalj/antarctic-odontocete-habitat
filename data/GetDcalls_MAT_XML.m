%% Compare Blue Whale D-call detections from MAT and XML files% This script:% 1. Loads detections from .mat files% 2. Loads detections from .xml files% 3. Concatenates detections separately% 4. Saves both to separate CSVs% 5. Computes comparison statistics% 6. Creates daily and weekly comparison plots

clear; clc;

%% ---------------- USER SETTINGS ----------------

site = 'SSI';

% MAT detector output 
foldermatDir = ['\snowman.ucsd.edu\Ally_Working_Disk\Analysis\Bm\Bm D call detector output\Antarc',site,'\data'];

% XML detector output
folderxmlDir = ['\snowman.ucsd.edu\Ally_Working_Disk\Analysis\Bm\Bm D call detector output\Antarc',site,'\submitted to tethys'];

% Output
folderoutDir = 'L:\Shared drives\Antarctic Marine Mammals\Marine Mammal Data\Mysticetes';

% Output 
CSVsmatCSV = fullfile(outDir, [site,'_BlueWhale_D_calls_MAT.csv']);xmlCSV = fullfile(outDir, [site,'_BlueWhale_D_calls_XML.csv']);

%% ---------------- LOAD MAT FILES ----------------

fprintf('\nLoading MAT detections...\n');

matFiles = dir(fullfile(matDir, '*.mat'));

matStart = datetime.empty(0,1);matEnd   = datetime.empty(0,1);matSrc   = strings(0,1);

for k = 1(matFiles)

fpath = fullfile(matFiles(k).folder, matFiles(k).name);

try
    load(fpath, '-mat');

    calls = hyd.detection.calls;

    [tStartRaw, tEndRaw] = getStartEnd(calls);

    tStart = toDatetime(tStartRaw);
    tEnd   = toDatetime(tEndRaw);

    n = min(numel(tStart), numel(tEnd));

    tStart = tStart(1:n);
    tEnd   = tEnd(1:n);

    matStart = [matStart; tStart]; %#ok<AGROW>
    matEnd   = [matEnd; tEnd]; %#ok<AGROW>
    matSrc   = [matSrc; repmat(string(matFiles(k).name), n, 1)]; %#ok<AGROW>

catch ME
    warning('Could not load MAT file: %s\n%s', matFiles(k).name, ME.message);
end

end

%% ---------------- SAVE MAT CSV ----------------

Tmat = table(matStart, matEnd, matSrc, ...'VariableNames', {'StartTime','EndTime','SourceFile'});

Tmat = sortrows(Tmat, 'StartTime');

writetable(Tmat, matCSV);

fprintf('Saved MAT detections: %d calls\n', height(Tmat));

%% ---------------- LOAD XML FILES ----------------

fprintf('\nLoading XML detections...\n');

xmlFiles = dir(fullfile(xmlDir, '*.xml'));

xmlStart = datetime.empty(0,1);xmlEnd   = datetime.empty(0,1);xmlSrc   = strings(0,1);

for k = 1(xmlFiles)

fpath = fullfile(xmlFiles(k).folder, xmlFiles(k).name);

try

    xDoc = xmlread(fpath);

    % ---- Find all Detection nodes ----
    detections = xDoc.getElementsByTagName('*');

    fileStarts = datetime.empty(0,1);
    fileEnds   = datetime.empty(0,1);

    for i = 0:detections.getLength-1

        node = detections.item(i);

        % Try common field names
        attrs = node.getAttributes;

        if isempty(attrs)
            continue;
        end

        startStr = '';
        endStr   = '';

        for j = 0:attrs.getLength-1

            attrName  = char(attrs.item(j).getName);
            attrValue = char(attrs.item(j).getValue);

            switch lower(attrName)

                case {'start','starttime','start_time'}
                    startStr = attrValue;

                case {'end','endtime','end_time'}
                    endStr = attrValue;
            end
        end

        % Convert if found
        if ~isempty(startStr)

            try
                t1 = parseXMLTime(startStr);

                if ~isempty(endStr)
                    t2 = parseXMLTime(endStr);
                else
                    t2 = t1;
                end

                fileStarts(end+1,1) = t1; %#ok<SAGROW>
                fileEnds(end+1,1)   = t2; %#ok<SAGROW>

            catch
            end
        end
    end

    n = min(numel(fileStarts), numel(fileEnds));

    xmlStart = [xmlStart; fileStarts(1:n)]; %#ok<AGROW>
    xmlEnd   = [xmlEnd; fileEnds(1:n)]; %#ok<AGROW>
    xmlSrc   = [xmlSrc; repmat(string(xmlFiles(k).name), n, 1)]; %#ok<AGROW>

catch ME
    warning('Could not load XML file: %s\n%s', xmlFiles(k).name, ME.message);
end

end

%% ---------------- SAVE XML CSV ----------------

Txml = table(xmlStart, xmlEnd, xmlSrc, ...'VariableNames', {'StartTime','EndTime','SourceFile'});

Txml = sortrows(Txml, 'StartTime');

writetable(Txml, xmlCSV);

fprintf('Saved XML detections: %d calls\n', height(Txml));

%% ---------------- BASIC COMPARISON STATS ----------------

fprintf('\n====================\n');fprintf('COMPARISON SUMMARY\n');fprintf('====================\n');

fprintf('MAT detections: %d\n', height(Tmat));fprintf('XML detections: %d\n', height(Txml));

fprintf('\nDate ranges:\n');

fprintf('MAT: %s --> %s\n', ...datestr(min(Tmat.StartTime)), ...datestr(max(Tmat.StartTime)));

fprintf('XML: %s --> %s\n', ...datestr(min(Txml.StartTime)), ...datestr(max(Txml.StartTime)));

%% ---------------- DAILY COUNTS ----------------

matDaily = groupsummary( ...table(dateshift(Tmat.StartTime,'start','day')), ...'Var1');

xmlDaily = groupsummary( ...table(dateshift(Txml.StartTime,'start','day')), ...'Var1');

matDaily.Properties.VariableNames = {'Date','GroupCount'};xmlDaily.Properties.VariableNames = {'Date','GroupCount'};

%% ---------------- WEEKLY COUNTS ----------------

matWeekly = groupsummary( ...table(dateshift(Tmat.StartTime,'start','week')), ...'Var1');

xmlWeekly = groupsummary( ...table(dateshift(Txml.StartTime,'start','week')), ...'Var1');

matWeekly.Properties.VariableNames = {'Date','GroupCount'};xmlWeekly.Properties.VariableNames = {'Date','GroupCount'};

%% ---------------- DAILY PLOT ----------------

figure('Color','w');

plot(matDaily.Date, matDaily.GroupCount, ...'LineWidth',1.5);

hold on

plot(xmlDaily.Date, xmlDaily.GroupCount, ...'LineWidth',1.5);

xlabel('Date');ylabel('Daily detections');

title([site ' Blue Whale D-calls: Daily Comparison']);

legend('MAT','XML');

grid on

saveas(gcf, fullfile(outDir, ...[site '_BlueWhale_DailyComparison.png']));

%% ---------------- WEEKLY PLOT ----------------

figure('Color','w');

plot(matWeekly.Date, matWeekly.GroupCount, ...'LineWidth',1.5);

hold on

plot(xmlWeekly.Date, xmlWeekly.GroupCount, ...'LineWidth',1.5);

xlabel('Date');ylabel('Weekly detections');

title([site ' Blue Whale D-calls: Weekly Comparison']);

legend('MAT','XML');

grid on

saveas(gcf, fullfile(outDir, ...[site '_BlueWhale_WeeklyComparison.png']));

fprintf('\nFinished!\n');

%% =========================================================%% ================= HELPER FUNCTIONS ======================%% =========================================================

function [tStartRaw, tEndRaw] = getStartEnd(calls)

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

error('Unsupported calls format.');

end

%% ---------------------------------------------------------

function dt = toDatetime(x)

if isempty(x)
    dt = datetime.empty;
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

        dt = datetime(x, ...
            'ConvertFrom','posixtime', ...
            'TimeZone','UTC');

    elseif xmax > 5e5

        dt = datetime(x, ...
            'ConvertFrom','datenum', ...
            'TimeZone','UTC');

    else

        error('Unknown numeric time format.');
    end

    return
end

error('Unsupported time type.');

end

%% ---------------------------------------------------------

function dt = parseXMLTime(tstr)

% Attempt ISO format first
try
    dt = datetime(tstr, ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss', ...
        'TimeZone','UTC');
    return
end

% Attempt fractional seconds
try
    dt = datetime(tstr, ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS', ...
        'TimeZone','UTC');
    return
end

% Generic fallback
dt = datetime(tstr, 'TimeZone','UTC');

end