%% Batch extract daily Fin/Blue 20 Hz detections from Tethys XML -> CSV
% Configure paths
inDir  = 'L:\Shared drives\Antarctic Marine Mammals\Marine Mammal Data\Mysticetes\xml_ignore\Bp_SSI01';
outDir = 'L:\Shared drives\Antarctic Marine Mammals\Marine Mammal Data\Mysticetes';
site   = 'SSI';
outCSV = fullfile(outDir, ['Antarc_',site,'_01_Bp_20Hz.csv']); % name as you like

% Scan files
files = dir(fullfile(inDir, '*.xml'));
allStart = datetime.empty(0,1);  allStart.TimeZone = 'UTC';
allEnd   = datetime.empty(0,1);  allEnd.TimeZone   = 'UTC';
allFile  = strings(0,1);
allScore = nan(0,1);
allCall  = strings(0,1);
allSpID  = strings(0,1);

for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);

    % --- Parse XML (prefer readstruct; fallback to DOM if needed) ---
    try
        X = readstruct(fpath, 'TextType','string');  % R2020b+
        % Expect fields: X.DataSource, X.Effort, X.OnEffort.Detection(:)
        dets = X.OnEffort.Detection;
        binMin = getBinMinutes(X);  % from Effort.Kind.Granularity @BinSize_m, default 1440
        [tStart, tEnd, score, call, spid] = extractDetections_struct(dets, binMin);
    catch
        % Fallback for older MATLAB or unexpected namespace quirks
        [tStart, tEnd, score, call, spid] = extractDetections_dom(fpath);
    end

    if isempty(tStart), continue; end

    allStart = [allStart; tStart]; %#ok<AGROW>
    allEnd   = [allEnd;   tEnd];   %#ok<AGROW>
    allFile  = [allFile;  repmat(string(files(k).name), numel(tStart), 1)]; %#ok<AGROW>
    allScore = [allScore; score]; %#ok<AGROW>
    allCall  = [allCall;  call];  %#ok<AGROW>
    allSpID  = [allSpID;  spid];  %#ok<AGROW>
end

% Assemble, sort, export
T = table(allStart, allScore, ...
    'VariableNames', {'StartTime','Score'});
T = sortrows(T, 'StartTime');
writetable(T, outCSV);
fprintf('Wrote %d rows to %s\n', height(T), outCSV);

%% ----------------- Helpers (keep at end) -----------------

function binMin = getBinMinutes(X)
% Read Effort.Kind.Granularity @BinSize_m; default to 1440 if missing
    binMin = 1440; % default: daily bins
    try
        G = X.Effort.Kind.Granularity;
        if isfield(G,'Attributes') && isfield(G.Attributes,'BinSize_m')
            binMin = str2double(G.Attributes.BinSize_m);
        elseif isfield(G,'BinSize_m')
            binMin = str2double(G.BinSize_m);
        end
        if isnan(binMin) || binMin<=0, binMin = 1440; end
    catch
        % keep default
    end
end

function [tStart, tEnd, score, call, spid] = extractDetections_struct(dets, binMin)
% Handle struct/struct-array as produced by readstruct
    if isempty(dets)
        tStart = datetime.empty(0,1); tStart.TimeZone = 'UTC';
        tEnd   = datetime.empty(0,1); tEnd.TimeZone   = 'UTC';
        score = nan(0,1); call = strings(0,1); spid = strings(0,1);
        return
    end

    if ~isstruct(dets)
        error('Unexpected OnEffort.Detection type: %s', class(dets));
    end

    % Normalize to an array of structs
    if ~isvector(dets), dets = dets(:); end
    n = numel(dets);

    tStart = NaT(n,1); tStart.TimeZone = 'UTC';
    tEnd   = NaT(n,1); tEnd.TimeZone   = 'UTC';
    score  = nan(n,1);
    call   = strings(n,1);
    spid   = strings(n,1);

    for i = 1:n
        % Start time as ISO 8601 (e.g., 2016-07-23T00:00:00.000Z)
        startStr = string(dets(i).Start);
        tStart(i) = iso2dt(startStr);
        tEnd(i)   = tStart(i) + minutes(binMin);

        % Optional fields
        try
            call(i) = string(dets(i).Call);
        catch, call(i) = ""; end
        try
            spid(i) = string(dets(i).SpeciesID);
        catch, spid(i) = ""; end

        try
            % Score nested under Parameters.Score
            if isfield(dets(i),'Parameters') && isfield(dets(i).Parameters,'Score')
                score(i) = double(dets(i).Parameters.Score);
            end
        catch
            % leave NaN
        end
    end
end

function [tStart, tEnd, score, call, spid] = extractDetections_dom(fpath)
% Namespace-tolerant DOM reader (older MATLAB). Computes End = Start + 1440 min if no bin size found.
    score = nan(0,1); call = strings(0,1); spid = strings(0,1);
    tStart = datetime.empty(0,1); tStart.TimeZone = 'UTC';
    tEnd   = datetime.empty(0,1); tEnd.TimeZone   = 'UTC';

    doc = xmlread(fpath);
    factory = javax.xml.xpath.XPathFactory.newInstance();
    xpath = factory.newXPath();

    % Helper to eval XPath and return a NodeList
    function nodes = xp(query)
        nodes = xpath.evaluate(query, doc, javax.xml.xpath.XPathConstants.NODESET);
    end

    % Bin size (minutes) if available; default 1440
    binMin = 1440;
    try
        attr = xpath.evaluate('string(//Granularity/@BinSize_m)', doc);
        if ~isempty(char(attr))
            v = str2double(char(attr));
            if ~isnan(v) && v>0, binMin = v; end
        end
    catch
    end

    detNodes = xp('//*[local-name()="OnEffort"]/*[local-name()="Detection"]');
    n = detNodes.getLength();
    if n==0, return; end

    tStart = NaT(n,1); tStart.TimeZone = 'UTC';
    tEnd   = NaT(n,1); tEnd.TimeZone   = 'UTC';
    score  = nan(n,1);
    call   = strings(n,1);
    spid   = strings(n,1);

    for i = 0:n-1
        det = detNodes.item(i);
        % Start
        s = char(xpath.evaluate('./*[local-name()="Start"]/text()', det));
        tStart(i+1) = iso2dt(string(s));
        tEnd(i+1)   = tStart(i+1) + minutes(binMin);
        % Optional fields
        c = char(xpath.evaluate('./*[local-name()="Call"]/text()', det));
        if ~isempty(c), call(i+1) = string(c); end
        id = char(xpath.evaluate('./*[local-name()="SpeciesID"]/text()', det));
        if ~isempty(id), spid(i+1) = string(id); end
        sc = char(xpath.evaluate('./*[local-name()="Parameters"]/*[local-name()="Score"]/text()', det));
        if ~isempty(sc), score(i+1) = str2double(sc); end
    end
end

function dt = iso2dt(s)
% Parse ISO 8601 like 2016-07-23T00:00:00.000Z into datetime UTC
    s = string(s);
    try
        % with milliseconds
        dt = datetime(s, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''', 'TimeZone','UTC');
    catch
        % without milliseconds
        dt = datetime(s, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''', 'TimeZone','UTC');
    end
end
