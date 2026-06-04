%% =========================================================
% Compare OLD verified blue whale D-call detections with NEW
% reprocessed detector output, estimate daily verification rate,
% and apply verification adjustment ONLY to missing/gap periods.
%
% Key logic:
% - Full verified/no-gap days:
%       FinalMetric = OldMetric
%
% - Partial gap days:
%       FinalMetric = OldMetric + adjusted NEW detections inside gap only
%
% - Completely unverified gap days:
%       FinalMetric = adjusted NEW detections inside gap only
%
% - Verification rate is estimated from nearest full verified day before
%   and nearest full verified day after each gap day.
%
% - If both nearest full verified days have OLD = 0, verification rate = 0.
%
% Final product is daily sums:
%       metricType = 'count'   -> number of detections per day
%       metricType = 'minutes' -> total detection minutes per day
%% =========================================================

clear;
clc;
close all;

%% =========================================================
%% USER SETTINGS
%% =========================================================

site = 'EI';

% Choose daily metric
% Options: 'count' or 'minutes'
metricType = 'count';

% Tolerance for matching OLD and NEW detections
matchTolerance = seconds(5);

switch site
    case 'SSI'
    % Date cutoff
    cutoff = datetime(2015,2,10);
    case 'EI'
    cutoff = datetime (2014,03,05);    
end

% Input files
oldCSV = ['L:\Shared drives\Antarctic Marine Mammals\Marine Mammal Data\Mysticetes\',site,'_BlueWhale_D_calls_XML.csv'];
newCSV = ['F:\Antarctic_Mysticetes\D calls\',site,'\filtered\All_Detections_Concatenated.csv'];

% Gap file
gapCSV = ['F:\Antarctic_Mysticetes\D calls\',site,'_gaps_only\',site,'_GapSummary.csv'];

% Output directory
outDir = ['F:\Antarctic_Mysticetes\D calls\',site,'\filtered'];

% Output files
dailyOutCSV = fullfile(outDir, [site '_Daily_Verification_Adjusted_Dcalls.csv']);
recoveredCSV = fullfile(outDir, [site '_Recovered_Dcalls_NeedValidation.csv']);

comparisonFig = fullfile(outDir, [site '_OldVsNew_DailyComparison.png']);
adjustedFig = fullfile(outDir, [site '_Daily_Verification_Adjusted_Dcalls.png']);
timelineFig = fullfile(outDir, [site '_OldVsNew_Timeline.png']);
nearestDiffFig = fullfile(outDir, [site '_NearestDetectionTimeDifferences.png']);

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% =========================================================
%% LOAD OLD AND NEW DETECTIONS
%% =========================================================

fprintf('\nLoading OLD detections...\n');
old = readtable(oldCSV);

fprintf('Loading NEW detections...\n');
new = readtable(newCSV);

%% =========================================================
%% STANDARDIZE DATETIME COLUMNS
%% =========================================================
old.StartTime = datetime(old.StartTime);
old.EndTime   = datetime(old.EndTime);

if any(strcmp(new.Properties.VariableNames, 'StartTime'))
    new.StartTime = datetime(new.StartTime);
else
    new.StartTime = datetime(new.Var1);
end

if any(strcmp(new.Properties.VariableNames, 'EndTime'))
    new.EndTime = datetime(new.EndTime);
else
    new.EndTime = datetime(new.Var2);
end

old.StartTime.TimeZone = '';
old.EndTime.TimeZone   = '';
new.StartTime.TimeZone = '';
new.EndTime.TimeZone   = '';

%% =========================================================
%% FILTER TO ANALYSIS PERIOD
%% =========================================================
old = old(old.StartTime >= cutoff, :);
new = new(new.StartTime >= cutoff, :);

fprintf('\n=== BASIC STATS AFTER CUTOFF ===\n');
fprintf('OLD detections: %d\n', height(old));
fprintf('NEW detections: %d\n', height(new));

fprintf('\nTime range:\n');
fprintf('OLD: %s to %s\n', datestr(min(old.StartTime)), datestr(max(old.StartTime)));
fprintf('NEW: %s to %s\n', datestr(min(new.StartTime)), datestr(max(new.StartTime)));

%% =========================================================
%% TIMELINE PLOT
%% =========================================================

figure('Color','w');
hold on;

plot(old.StartTime, ones(height(old),1), '.r');
plot(new.StartTime, 2*ones(height(new),1), '.k');

yticks([1 2]);
yticklabels({'OLD verified','NEW reprocessed'});
xlabel('Date');
title([site ' Blue Whale D-calls: Detection Timeline Comparison']);
grid on;

saveas(gcf, timelineFig);

%% =========================================================
%% DAILY COUNTS / MINUTES FOR OLD AND NEW
%% =========================================================

old.Date = dateshift(old.StartTime, 'start', 'day');
new.Date = dateshift(new.StartTime, 'start', 'day');

switch lower(metricType)

    case 'count'
        old.Metric = ones(height(old),1);
        new.Metric = ones(height(new),1);
        ylab = 'Daily detections';

    case 'minutes'
        old.Metric = minutes(old.EndTime - old.StartTime);
        new.Metric = minutes(new.EndTime - new.StartTime);
        ylab = 'Daily detection minutes';

    otherwise
        error('metricType must be either count or minutes.');
end

[G_old, oldDateGroup] = findgroups(old.Date);
oldMetricSum = splitapply(@sum, old.Metric, G_old);
oldDaily = table(oldDateGroup, oldMetricSum, ...
    'VariableNames', {'Date','OldMetric'});

[G_new, newDateGroup] = findgroups(new.Date);
newMetricSum = splitapply(@sum, new.Metric, G_new);
newDaily = table(newDateGroup, newMetricSum, ...
    'VariableNames', {'Date','NewMetric'});

allDates = unique([oldDaily.Date; newDaily.Date]);
Daily = table(allDates, 'VariableNames', {'Date'});

Daily = outerjoin(Daily, oldDaily, ...
    'Keys','Date', ...
    'MergeKeys',true);

Daily = outerjoin(Daily, newDaily, ...
    'Keys','Date', ...
    'MergeKeys',true);

Daily.OldMetric(isnan(Daily.OldMetric)) = 0;
Daily.NewMetric(isnan(Daily.NewMetric)) = 0;

%% =========================================================
%% BASIC DAILY OLD VS NEW PLOT
%% =========================================================

figure('Color','w');
hold on;

plot(Daily.Date, Daily.OldMetric, '-r', 'LineWidth', 1.5);
plot(Daily.Date, Daily.NewMetric, '-k', 'LineWidth', 1.5);

legend('OLD verified','NEW reprocessed');
xlabel('Date');
ylabel(ylab);
title([site ' Blue Whale D-calls: Daily OLD vs NEW']);
grid on;

saveas(gcf, comparisonFig);

%% =========================================================
%% MATCH NEW DETECTIONS TO OLD VERIFIED DETECTIONS
%% =========================================================

fprintf('\nMatching NEW detections to OLD verified detections...\n');

isMatchedNew = false(height(new),1);

for i = 1:height(new)

    t = new.StartTime(i);

    isMatchedNew(i) = any(abs(old.StartTime - t) < matchTolerance);
end

new.VerificationStatus = strings(height(new),1);
new.VerificationStatus(isMatchedNew)  = "Previously verified";
new.VerificationStatus(~isMatchedNew) = "Recovered / unverified";

verified_new = new(isMatchedNew,:);
recovered = new(~isMatchedNew,:);

fprintf('\n=== RECOVERED / UNVERIFIED DETECTIONS ===\n');
fprintf('Previously verified detections in NEW: %d\n', height(verified_new));
fprintf('Recovered / unverified detections in NEW: %d\n', height(recovered));
fprintf('Percent recovered / unverified: %.2f%%\n', ...
    100 * height(recovered) / height(new));

writetable(recovered, recoveredCSV);

fprintf('\nSaved recovered detections:\n%s\n', recoveredCSV);

%% =========================================================
%% NEAREST DETECTION TIME DIFFERENCES
%% =========================================================

fprintf('\nComputing nearest OLD-to-NEW detection time differences...\n');

diffs = [];

for i = 1:height(old)

    t = old.StartTime(i);

    if height(new) > 0
        d = min(abs(new.StartTime - t));
        diffs(end+1,1) = seconds(d); %#ok<SAGROW>
    end
end

figure('Color','w');
histogram(diffs, 50);
xlabel('Nearest OLD-to-NEW detection time difference (seconds)');
ylabel('Count');
title([site ' Nearest Detection Time Differences']);
grid on;

saveas(gcf, nearestDiffFig);

matchedOld = sum(diffs < seconds(matchTolerance));

fprintf('\nMatched OLD detections within %.1f sec: %d / %d\n', ...
    seconds(matchTolerance), matchedOld, height(old));

%% =========================================================
%% LOAD GAP TABLE AND CLASSIFY DAYS
%% =========================================================

fprintf('\nLoading gap table...\n');

gaps = readtable(gapCSV);

gaps.GapStart = datetime(gaps.GapStart);
gaps.GapEnd   = datetime(gaps.GapEnd);

gaps.GapStart.TimeZone = '';
gaps.GapEnd.TimeZone   = '';

Daily.GapStatus = strings(height(Daily),1);
Daily.GapStatus(:) = "Full verified day";

Daily.GapMinutes = zeros(height(Daily),1);

for i = 1:height(Daily)

    dayStart = Daily.Date(i);
    dayEnd   = dayStart + days(1);

    gapOverlapIdx = find(gaps.GapStart < dayEnd & gaps.GapEnd > dayStart);

    totalGapMinutes = 0;

    for g = gapOverlapIdx'

        thisGapStart = max(gaps.GapStart(g), dayStart);
        thisGapEnd   = min(gaps.GapEnd(g), dayEnd);

        totalGapMinutes = totalGapMinutes + minutes(thisGapEnd - thisGapStart);
    end

    Daily.GapMinutes(i) = totalGapMinutes;

    if totalGapMinutes >= 1440
        Daily.GapStatus(i) = "Completely unverified gap day";
    elseif totalGapMinutes > 0
        Daily.GapStatus(i) = "Partial gap day";
    end
end

%% =========================================================
%% ESTIMATE VERIFICATION RATE FROM FULL VERIFIED DAYS
%% =========================================================

Daily.VerificationRate = NaN(height(Daily),1);

fullIdx = Daily.GapStatus == "Full verified day" & Daily.NewMetric > 0;

Daily.VerificationRate(fullIdx) = Daily.OldMetric(fullIdx) ./ Daily.NewMetric(fullIdx);

Daily.VerificationRate(Daily.VerificationRate > 1) = 1;

globalVerificationRate = median(Daily.VerificationRate(fullIdx), 'omitnan');

fprintf('\n=== VERIFICATION RATE ESTIMATION ===\n');
fprintf('Number of full verified days used: %d\n', sum(fullIdx));
fprintf('Median full-day verification rate: %.2f%%\n', ...
    100 * globalVerificationRate);

%% =========================================================
%% ASSIGN LOCAL VERIFICATION RATE TO GAP DAYS
%% =========================================================

gapIdx = Daily.GapStatus ~= "Full verified day";

for i = find(gapIdx)'

    beforeIdx = find(fullIdx & Daily.Date < Daily.Date(i), 1, 'last');
    afterIdx  = find(fullIdx & Daily.Date > Daily.Date(i), 1, 'first');

    neighborRates = [];
    neighborOldMetrics = [];

    if ~isempty(beforeIdx)
        neighborRates(end+1) = Daily.VerificationRate(beforeIdx); %#ok<SAGROW>
        neighborOldMetrics(end+1) = Daily.OldMetric(beforeIdx); %#ok<SAGROW>
    end

    if ~isempty(afterIdx)
        neighborRates(end+1) = Daily.VerificationRate(afterIdx); %#ok<SAGROW>
        neighborOldMetrics(end+1) = Daily.OldMetric(afterIdx); %#ok<SAGROW>
    end

    if ~isempty(neighborOldMetrics) && all(neighborOldMetrics == 0)
        Daily.VerificationRate(i) = 0;

    elseif ~isempty(neighborRates)
        Daily.VerificationRate(i) = mean(neighborRates, 'omitnan');

    else
        Daily.VerificationRate(i) = globalVerificationRate;
    end
end

Daily.VerificationRate(isnan(Daily.VerificationRate)) = globalVerificationRate;

%% =========================================================
%% APPLY VERIFICATION RATE ONLY TO NEW DETECTIONS INSIDE GAPS
%% =========================================================

Daily.NewInGapMetric = zeros(height(Daily),1);
Daily.VerificationAdjustedGapMetric = zeros(height(Daily),1);

% Default: use old verified metric
Daily.FinalMetric = Daily.OldMetric;

for i = 1:height(Daily)

    dayStart = Daily.Date(i);
    dayEnd   = dayStart + days(1);

    if Daily.GapStatus(i) == "Full verified day"
        continue
    end

    gapOverlapIdx = find(gaps.GapStart < dayEnd & gaps.GapEnd > dayStart);

    newInGapMetric = 0;

    for g = gapOverlapIdx'

        thisGapStart = max(gaps.GapStart(g), dayStart);
        thisGapEnd   = min(gaps.GapEnd(g), dayEnd);

        inGap = new.StartTime >= thisGapStart & new.StartTime < thisGapEnd;

        switch lower(metricType)

            case 'count'
                newInGapMetric = newInGapMetric + sum(inGap);

            case 'minutes'
                newInGapMetric = newInGapMetric + ...
                    sum(minutes(new.EndTime(inGap) - new.StartTime(inGap)));
        end
    end

    Daily.NewInGapMetric(i) = newInGapMetric;

    Daily.VerificationAdjustedGapMetric(i) = ...
        Daily.NewInGapMetric(i) * Daily.VerificationRate(i);

    % Partial gap days:
    % Final = old verified detections outside gap + adjusted new detections inside gap
    %
    % Completely unverified gap days:
    % OldMetric should usually be zero, so this becomes adjusted gap detections
    Daily.FinalMetric(i) = ...
        Daily.OldMetric(i) + Daily.VerificationAdjustedGapMetric(i);
end

if strcmpi(metricType, 'count')
    Daily.NewInGapMetric = round(Daily.NewInGapMetric);
    Daily.VerificationAdjustedGapMetric = round(Daily.VerificationAdjustedGapMetric);
    Daily.FinalMetric = round(Daily.FinalMetric);
end

%% =========================================================
%% SAVE FINAL DAILY TABLE
%% =========================================================

writetable(Daily, dailyOutCSV);

fprintf('\nSaved daily verification-adjusted table:\n%s\n', dailyOutCSV);

%% =========================================================
%% PLOT VERIFICATION-ADJUSTED DAILY VALUES
%% =========================================================

figure('Color','w');
hold on;

plot(Daily.Date, Daily.OldMetric, '-r', 'LineWidth', 1.5);
plot(Daily.Date, Daily.NewMetric, '-k', 'LineWidth', 1.5);
plot(Daily.Date, Daily.FinalMetric, '-g', 'LineWidth', 1.5);

legend('OLD verified','NEW raw','NEW all-day verification-adjusted','Final daily metric');
xlabel('Date');
ylabel(ylab);
title([site ' Blue Whale D-calls: Daily Verification Adjustment']);
grid on;

saveas(gcf, adjustedFig);

%% =========================================================
%% PRINT SUMMARY BY GAP STATUS
%% =========================================================

fprintf('\n=== SUMMARY BY GAP STATUS ===\n');

statuses = unique(Daily.GapStatus);

for i = 1:numel(statuses)

    idx = Daily.GapStatus == statuses(i);

    fprintf('\n%s\n', statuses(i));
    fprintf('Days: %d\n', sum(idx));
    fprintf('Mean gap minutes: %.2f\n', mean(Daily.GapMinutes(idx), 'omitnan'));
    fprintf('Mean OLD metric: %.2f\n', mean(Daily.OldMetric(idx), 'omitnan'));
    fprintf('Mean NEW metric: %.2f\n', mean(Daily.NewMetric(idx), 'omitnan'));
    fprintf('Mean NEW in-gap metric: %.2f\n', mean(Daily.NewInGapMetric(idx), 'omitnan'));
    fprintf('Mean adjusted gap metric: %.2f\n', mean(Daily.VerificationAdjustedGapMetric(idx), 'omitnan'));
    fprintf('Mean final metric: %.2f\n', mean(Daily.FinalMetric(idx), 'omitnan'));
    fprintf('Mean verification rate: %.2f%%\n', ...
        100 * mean(Daily.VerificationRate(idx), 'omitnan'));
end

fprintf('\nFinished!\n');