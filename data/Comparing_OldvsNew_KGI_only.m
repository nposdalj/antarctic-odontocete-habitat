old = readtable('L:\Shared drives\Antarctic Marine Mammals\Marine Mammal Data\Mysticetes\EI_BlueWhale_D_calls_XML.csv');
new = readtable('F:\Antarctic_Mysticetes\D calls\EI\filtered\All_Detections_Concatenated.csv');

% Ensure datetime
old.StartTime = datetime(old.StartTime);
old.EndTime   = datetime(old.EndTime);

new.StartTime = datetime(new.Var1);
new.EndTime   = datetime(new.Var2);

% ==============================
% FILTER: Only keep >= Feb 10, 2015
% ==============================
cutoff = datetime(2015,2,10);

old = old(old.StartTime >= cutoff, :);
new = new(new.StartTime >= cutoff, :);

fprintf('\n=== BASIC STATS ===\n');

fprintf('OLD detections: %d\n', height(old));
fprintf('NEW detections: %d\n', height(new));

fprintf('\nTime range:\n');
fprintf('OLD: %s to %s\n', datestr(min(old.StartTime)), datestr(max(old.StartTime)));
fprintf('NEW: %s to %s\n', datestr(min(new.StartTime)), datestr(max(new.StartTime)));

figure; hold on

plot(old.StartTime, ones(height(old),1), '.r')
plot(new.StartTime, 2*ones(height(new),1), '.k')

yticks([1 2])
yticklabels({'OLD','NEW'})
xlabel('Date')
title('Detection Timeline Comparison')
grid on

% Bin to daily counts
old_day = dateshift(old.StartTime, 'start', 'day');
new_day = dateshift(new.StartTime, 'start', 'day');

old_counts = groupsummary(table(old_day), 'old_day');
new_counts = groupsummary(table(new_day), 'new_day');

figure; hold on
plot(old_counts.old_day, old_counts.GroupCount, '-r', 'LineWidth', 1.5)
plot(new_counts.new_day, new_counts.GroupCount, '-k', 'LineWidth', 1.5)

legend('OLD','NEW')
xlabel('Date')
ylabel('Detections per day')
title('Daily Detection Counts')
grid on

fprintf('\n=== GAP ANALYSIS ===\n');

old_dt = diff(sort(old.StartTime));
new_dt = diff(sort(new.StartTime));

fprintf('Largest OLD gap: %.2f days\n', max(days(old_dt)));
fprintf('Largest NEW gap: %.2f days\n', max(days(new_dt)));

tol = seconds(5); % tolerance window

matches = 0;

for i = 1:height(old)
    
    t = old.StartTime(i);
    
    if any(abs(new.StartTime - t) < tol)
        matches = matches + 1;
    end
end

fprintf('\nMatched detections within 5 sec: %d / %d\n', matches, height(old));

diffs = [];

for i = 1:height(old)
    
    t = old.StartTime(i);
    [d, idx] = min(abs(new.StartTime - t));
    
    diffs(end+1) = seconds(d);
end

figure;
histogram(diffs, 50)
xlabel('Time difference (seconds)')
title('Nearest Detection Time Differences')