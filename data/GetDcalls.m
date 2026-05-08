%% Compare Blue Whale D-call detections from MAT and XML files

clear; clc;

%% --------------------------------------------------------
%% USER SETTINGS
%% --------------------------------------------------------

site = 'SSI';

% MAT folder
matDir = ['\\snowman.ucsd.edu\Ally_Working_Disk\Analysis\Bm\Bm D call detector output\Antarc\',site,'\Data'];

% XML folder
xmlDir = ['\\snowman.ucsd.edu\Ally_Working_Disk\Analysis\Bm\Bm D call detector output\Antarc\',site,'\submitted to tethys'];

% Output folder
outDir = 'L:\Shared drives\Antarctic Marine Mammals\Marine Mammal Data\Mysticetes';

if ~exist(outDir,'dir')
    mkdir(outDir)
end

matCSV = fullfile(outDir,[site '_BlueWhale_D_calls_MAT.csv']);
xmlCSV = fullfile(outDir,[site '_BlueWhale_D_calls_XML.csv']);

%% --------------------------------------------------------
%% LOAD MAT FILES
%% --------------------------------------------------------

fprintf('\nLoading MAT detections...\n')

matFiles = dir(fullfile(matDir,'*.mat'));

allStartMAT = datetime.empty(0,1);
allStartMAT.TimeZone = 'UTC';

allEndMAT = datetime.empty(0,1);
allEndMAT.TimeZone = 'UTC';

allSrcMAT = strings(0,1);

for k = 1:numel(matFiles)

    try

        fpath = fullfile(matFiles(k).folder, matFiles(k).name);

        load(fpath,'-mat');

        calls = hyd.detection.calls;

        [tStartRaw,tEndRaw] = getStartEnd(calls);

        tStart = toDatetime(tStartRaw);
        tEnd   = toDatetime(tEndRaw);

        n = min(numel(tStart), numel(tEnd));

        tStart = tStart(1:n);
        tEnd   = tEnd(1:n);

        allStartMAT = [allStartMAT; tStart];
        allEndMAT   = [allEndMAT; tEnd];

        allSrcMAT = [allSrcMAT; ...
            repmat(string(matFiles(k).name),n,1)];

    catch ME

        warning('Problem reading MAT file: %s\n%s', ...
            matFiles(k).name, ME.message);

    end
end

%% Save MAT CSV

Tmat = table( ...
    allStartMAT, ...
    allEndMAT, ...
    allSrcMAT, ...
    'VariableNames',{'StartTime','EndTime','SourceFile'});

Tmat = sortrows(Tmat,'StartTime');

writetable(Tmat, matCSV);

fprintf('Saved MAT detections: %d\n', height(Tmat));

%% --------------------------------------------------------
%% LOAD XML FILES
%% --------------------------------------------------------

fprintf('\nLoading XML detections...\n')

xmlFiles = dir(fullfile(xmlDir,'*.xml'));

allStartXML = datetime.empty(0,1);
allStartXML.TimeZone = 'UTC';

allEndXML = datetime.empty(0,1);
allEndXML.TimeZone = 'UTC';

allSrcXML = strings(0,1);

for k = 1:numel(xmlFiles)

    try

        fpath = fullfile(xmlFiles(k).folder, xmlFiles(k).name);

        xDoc = xmlread(fpath);

        detections = xDoc.getElementsByTagName('Detection');

        nDet = detections.getLength;

        for i = 0:nDet-1

            det = detections.item(i);

            children = det.getChildNodes;

            startStr = '';
            endStr   = '';

            for j = 0:children.getLength-1

                node = children.item(j);

                if node.getNodeType ~= node.ELEMENT_NODE
                    continue
                end

                nodeName = char(node.getNodeName);

                switch nodeName

                    case 'Start'
                        startStr = char(node.getTextContent);

                    case 'End'
                        endStr = char(node.getTextContent);
                end
            end

            if ~isempty(startStr)

                tStart = datetime( ...
                    startStr, ...
                    'InputFormat','yyyy-MM-dd''T''HH:mm:ssX', ...
                    'TimeZone','UTC');

                tEnd = datetime( ...
                    endStr, ...
                    'InputFormat','yyyy-MM-dd''T''HH:mm:ssX', ...
                    'TimeZone','UTC');

                allStartXML(end+1,1) = tStart;
                allEndXML(end+1,1)   = tEnd;

                allSrcXML(end+1,1) = string(xmlFiles(k).name);

            end
        end

    catch ME

        warning('Problem reading XML file: %s\n%s', ...
            xmlFiles(k).name, ME.message);

    end
end

%% Save XML CSV

Txml = table( ...
    allStartXML, ...
    allEndXML, ...
    allSrcXML, ...
    'VariableNames',{'StartTime','EndTime','SourceFile'});

Txml = sortrows(Txml,'StartTime');

writetable(Txml, xmlCSV);

fprintf('Saved XML detections: %d\n', height(Txml));

%% --------------------------------------------------------
%% SUMMARY STATS
%% --------------------------------------------------------

fprintf('\n====================\n')
fprintf('COMPARISON SUMMARY\n')
fprintf('====================\n')

fprintf('MAT detections: %d\n', height(Tmat));
fprintf('XML detections: %d\n', height(Txml));

if ~isempty(Tmat)
    fprintf('MAT range: %s --> %s\n', ...
        datestr(min(Tmat.StartTime)), ...
        datestr(max(Tmat.StartTime)));
end

if ~isempty(Txml)
    fprintf('XML range: %s --> %s\n', ...
        datestr(min(Txml.StartTime)), ...
        datestr(max(Txml.StartTime)));
end

%% --------------------------------------------------------
%% DAILY COUNTS
%% --------------------------------------------------------

if ~isempty(Tmat)

    matDaily = groupsummary( ...
        table(dateshift(Tmat.StartTime,'start','day')), ...
        'Var1');

    matDaily.Properties.VariableNames = ...
        {'Date','Count'};

else

    matDaily = table;

end

if ~isempty(Txml)

    xmlDaily = groupsummary( ...
        table(dateshift(Txml.StartTime,'start','day')), ...
        'Var1');

    xmlDaily.Properties.VariableNames = ...
        {'Date','Count'};

else

    xmlDaily = table;

end

%% --------------------------------------------------------
%% DAILY PLOT
%% --------------------------------------------------------

figure('Color','w');
hold on

if ~isempty(matDaily)

    plot(matDaily.Date, matDaily.Count, ...
        'LineWidth',1.5);

end

if ~isempty(xmlDaily)

    plot(xmlDaily.Date, xmlDaily.Count, ...
        'LineWidth',1.5);

end

xlabel('Date')
ylabel('Daily detections')

title([site ' Blue Whale D-call Daily Comparison'])

legendEntries = {};

if ~isempty(matDaily)
    legendEntries{end+1} = 'MAT';
end

if ~isempty(xmlDaily)
    legendEntries{end+1} = 'XML';
end

if ~isempty(legendEntries)
    legend(legendEntries)
end

grid on

saveas(gcf, fullfile(outDir, ...
    [site '_BlueWhale_DailyComparison.png']));

%% --------------------------------------------------------
%% WEEKLY COUNTS
%% --------------------------------------------------------

if ~isempty(Tmat)

    matWeekly = groupsummary( ...
        table(dateshift(Tmat.StartTime,'start','week')), ...
        'Var1');

    matWeekly.Properties.VariableNames = ...
        {'Date','Count'};

else

    matWeekly = table;

end

if ~isempty(Txml)

    xmlWeekly = groupsummary( ...
        table(dateshift(Txml.StartTime,'start','week')), ...
        'Var1');

    xmlWeekly.Properties.VariableNames = ...
        {'Date','Count'};

else

    xmlWeekly = table;

end

%% --------------------------------------------------------
%% WEEKLY PLOT
%% --------------------------------------------------------

figure('Color','w');
hold on

if ~isempty(matWeekly)

    plot(matWeekly.Date, matWeekly.Count, ...
        'LineWidth',1.5);

end

if ~isempty(xmlWeekly)

    plot(xmlWeekly.Date, xmlWeekly.Count, ...
        'LineWidth',1.5);

end

xlabel('Date')
ylabel('Weekly detections')

title([site ' Blue Whale D-call Weekly Comparison'])

legendEntries = {};

if ~isempty(matWeekly)
    legendEntries{end+1} = 'MAT';
end

if ~isempty(xmlWeekly)
    legendEntries{end+1} = 'XML';
end

if ~isempty(legendEntries)
    legend(legendEntries)
end

grid on

saveas(gcf, fullfile(outDir, ...
    [site '_BlueWhale_WeeklyComparison.png']));

fprintf('\nFinished!\n')

%% ========================================================
%% HELPER FUNCTIONS
%% ========================================================

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

    error('Unsupported calls type')
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

            error('Unknown numeric time format')

        end

        return
    end

    error('Unsupported time type')
end