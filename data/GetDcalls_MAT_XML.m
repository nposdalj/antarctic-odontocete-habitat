%% Compare Blue Whale D-call detections from MAT and XML files% This script:% 1. Loads detections from .mat files% 2. Loads detections from .xml files% 3. Concatenates detections separately% 4. Saves both to separate CSVs% 5. Computes comparison statistics% 6. Creates daily and weekly comparison plots

clear; clc;

%% ---------------- USER SETTINGS ----------------

site = 'EI';

% MAT detector output 
foldermatDir = ['\\snowman.ucsd.edu\Ally_Working_Disk\Analysis\Bm\Bm D call detector output\Antarc\',site,'\data'];

% XML detector output
folderxmlDir = ['\\snowman.ucsd.edu\Ally_Working_Disk\Analysis\Bm\Bm D call detector output\Antarc\',site,'\submitted to tethys'];

% Output
folderoutDir = 'L:\Shared drives\Antarctic Marine Mammals\Marine Mammal Data\Mysticetes';

% Output 
CSVsmatCSV = fullfile(folderoutDir, [site,'_BlueWhale_D_calls_MAT.csv']);
xmlCSV = fullfile(folderoutDir, [site,'_BlueWhale_D_calls_XML.csv']);

%% ---------------- LOAD MAT FILES ----------------

fprintf('\nLoading MAT detections...\n');

matFiles = dir(fullfile(foldermatDir, '*.mat'));

matStart = datetime.empty(0,1);
matEnd   = datetime.empty(0,1);matSrc   = strings(0,1);

for k = 1:length(matFiles)

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

    % Remove timezone inconsistencies
tStart.TimeZone = '';
tEnd.TimeZone   = '';

if ~isempty(matStart)
    matStart.TimeZone = '';
end

if ~isempty(matEnd)
    matEnd.TimeZone = '';
end

    matStart = [matStart; tStart]; %#ok<AGROW>
    matEnd   = [matEnd; tEnd]; %#ok<AGROW>
    matSrc   = [matSrc; repmat(string(matFiles(k).name), n, 1)]; %#ok<AGROW>

catch ME
    warning('Could not load MAT file: %s\n%s', matFiles(k).name, ME.message);
end

end

%% ---------------- SAVE MAT CSV ----------------

Tmat = table(matStart, matEnd, matSrc, ...
    'VariableNames', {'StartTime','EndTime','SourceFile'});

Tmat = sortrows(Tmat, 'StartTime');

writetable(Tmat, CSVsmatCSV);

fprintf('Saved MAT detections: %d calls\n', height(Tmat));

%% ---------------- LOAD XML FILES ----------------

fprintf('\nLoading XML detections...\n');

xmlFiles = dir(fullfile(folderxmlDir, '*.xml'));

xmlStart = datetime.empty(0,1);xmlEnd   = datetime.empty(0,1);xmlSrc   = strings(0,1);

for k = 1:length(xmlFiles)

fpath = fullfile(xmlFiles(k).folder, xmlFiles(k).name);

try

  xDoc = xmlread(fpath);

% Get all Start and End nodes
startNodes = xDoc.getElementsByTagName('Start');
endNodes   = xDoc.getElementsByTagName('End');

fileStarts = datetime.empty(0,1);
fileEnds   = datetime.empty(0,1);

nNodes = min(startNodes.getLength, endNodes.getLength);

for i = 0:nNodes-1

    try

        startStr = char(startNodes.item(i).getTextContent);
        endStr   = char(endNodes.item(i).getTextContent);

        t1 = parseXMLTime(startStr);
        t2 = parseXMLTime(endStr);

        % Remove timezone for consistency
        t1.TimeZone = '';
        t2.TimeZone = '';

        fileStarts(end+1,1) = t1; %#ok<SAGROW>
        fileEnds(end+1,1)   = t2; %#ok<SAGROW>

    catch ME

        warning('Problem parsing XML detection: %s', ME.message);

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

Txml = table(xmlStart, xmlEnd, xmlSrc, ...
    'VariableNames', {'StartTime','EndTime','SourceFile'});

Txml = sortrows(Txml, 'StartTime');

writetable(Txml, xmlCSV);

fprintf('Saved XML detections: %d calls\n', height(Txml));

%% ---------------- BASIC COMPARISON STATS ----------------

fprintf('\n====================\n');fprintf('COMPARISON SUMMARY\n');fprintf('====================\n');

fprintf('MAT detections: %d\n', height(Tmat));fprintf('XML detections: %d\n', height(Txml));

fprintf('\nDate ranges:\n');

fprintf('MAT: %s --> %s\n', ...
    datestr(min(Tmat.StartTime)), ...
    datestr(max(Tmat.StartTime)));

fprintf('XML: %s --> %s\n', ...
    datestr(min(Txml.StartTime)), ...
    datestr(max(Txml.StartTime)));

%% ---------------- DAILY COUNTS ----------------

matDaily = groupsummary( ...
    table(dateshift(Tmat.StartTime,'start','day')), ...
    'Var1');

xmlDaily = groupsummary( ...
    table(dateshift(Txml.StartTime,'start','day')), ...
    'Var1');

matDaily.Properties.VariableNames = {'Date','GroupCount'};
xmlDaily.Properties.VariableNames = {'Date','GroupCount'};

%% ---------------- WEEKLY COUNTS ----------------

matWeekly = groupsummary( ...
    table(dateshift(Tmat.StartTime,'start','week')), ...
    'Var1');

xmlWeekly = groupsummary( ...
    table(dateshift(Txml.StartTime,'start','week')), ...
    'Var1');

matWeekly.Properties.VariableNames = {'Date','GroupCount'};
xmlWeekly.Properties.VariableNames = {'Date','GroupCount'};

%% ---------------- DAILY PLOT ----------------

figure('Color','w');

plot(matDaily.Date, matDaily.GroupCount, ...
    'LineWidth',1.5);

hold on

plot(xmlDaily.Date, xmlDaily.GroupCount, ...
    'LineWidth',1.5);

xlabel('Date');ylabel('Daily detections');

title([site ' Blue Whale D-calls: Daily Comparison']);

legend('MAT','XML');

grid on

saveas(gcf, fullfile(folderoutDir, ...
    [site '_BlueWhale_DailyComparison.png']));

%% ---------------- WEEKLY PLOT ----------------

figure('Color','w');

plot(matWeekly.Date, matWeekly.GroupCount, ...
    'LineWidth',1.5);

hold on

plot(xmlWeekly.Date, xmlWeekly.GroupCount, ...
    'LineWidth',1.5);

xlabel('Date');ylabel('Weekly detections');

title([site ' Blue Whale D-calls: Weekly Comparison']);

legend('MAT','XML');

grid on

saveas(gcf, fullfile(folderoutDir, ...
    [site '_BlueWhale_WeeklyComparison.png']));

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

% Remove trailing Z
tstr = erase(tstr, 'Z');

try

    dt = datetime(tstr, ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS', ...
        'TimeZone','UTC');

catch

    dt = datetime(tstr, ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss', ...
        'TimeZone','UTC');

end

end