%% Script to convert daily .mat files to required format for processing

clearvars;

% Define input and output folders
inputFolder = 'C:\Users\nposd\Documents\GitHub\antarctic-odontocete-habitat\Environmental Data\HYCOM\Data'; % Your existing daily .mat files
outputFolder = 'C:\Users\nposd\Documents\GitHub\antarctic-odontocete-habitat\Environmental Data\HYCOM\Formatted_Daily_Files'; % Output formatted files

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% List daily .mat files in inputFolder
fileList = dir(fullfile(inputFolder, '*.mat'));

% Loop through each file and reformat
for k = 1:length(fileList)
    filename = fullfile(inputFolder, fileList(k).name);
    fprintf('Processing %s\n', filename);

    data = load(filename);

    % Define variables from loaded data
    lon = data.lon;
    lat = data.lat;
    depth = data.depth;
    time = data.time;

    % Extract variables (already daily averaged)
    dayMeanTemperature = data.temp;
    dayMeanSalinity = data.salt;
    dayMeanUVel = data.u;
    dayMeanVVel = data.v;

    % Check dimensions are correct (lon x lat x depth)
    % Adjust dimensions if necessary using permute/reshape

    % Handle missing values if needed (assuming 1000 as missing value)
    missValue = 1000;
    dayMeanTemperature(abs(dayMeanTemperature) >= missValue) = nan;
    dayMeanSalinity(abs(dayMeanSalinity) >= missValue) = nan;
    dayMeanUVel(abs(dayMeanUVel) >= missValue) = nan;
    dayMeanVVel(abs(dayMeanVVel) >= missValue) = nan;

    % Create Houridx as empty to match structure if required
    Houridx = {};

    % Save formatted data to new .mat file
    outputFilename = fullfile(outputFolder, [datestr(time,'yyyymmdd') '_daily.mat']);
    save(outputFilename, 'dayMeanTemperature', 'dayMeanSalinity', 'dayMeanUVel',...
        'dayMeanVVel', 'lon', 'lat', 'depth', 'time', 'Houridx');

end

fprintf('All files processed and saved to %s\n', outputFolder);
