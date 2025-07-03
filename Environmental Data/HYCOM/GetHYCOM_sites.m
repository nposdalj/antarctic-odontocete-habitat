%% Elephant Island
%{
clc;clearvars
aimpath = "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/HYCOM/EI";
region = [303.406 304.406 -61.3869 -59.3869]; % +/- 0.5 degrees from EI
timeList = datetime(2014,3,5):datetime(2014,7,17);
varList = {'ssh','temp','salt','u','v'};
URL = 'http://tds.hycom.org/thredds/catalogs/GLBv0.08/expt_53.X.html';

nTimes = numel(timeList);
for iTime = 1:nTimes
     timeTick = timeList(iTime);
     D = get_hycom_online(aimpath, region, timeTick, varList);
end

clc;clearvars
aimpath = "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/HYCOM/EI";
region = [303.406 304.406 -61.3869 -59.3869]; % +/- 0.5 degrees from EI
timeList = datetime(2014,3,5):datetime(2014,7,17);
varList = {'ssh','temp','salt','u','v'};
URL = 'http://tds.hycom.org/thredds/catalogs/GLBv0.08/expt_56.3.html';

nTimes = numel(timeList);
for iTime = 1:nTimes
     timeTick = timeList(iTime);
     D = get_hycom_online(aimpath, region, timeTick, varList);
end
%}

%% King George Island
clc;clearvars
aimpath = "C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/HYCOM/KGI";
region = [302.558083 303.558083 -61.957817 -60.957817]; % +/- 0.5 degrees from EI
timeList = datetime(2015,2,10):datetime(2016,1,29);
varList = {'ssh','temp','salt','u','v'};
URL = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_56.3?';

nTimes = numel(timeList);
for iTime = 1:nTimes
     timeTick = timeList(iTime);
     which get_hycom_online
     D = get_hycom_online(aimpath, region, timeTick, varList, URL);
end
