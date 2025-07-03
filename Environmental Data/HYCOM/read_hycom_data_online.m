clc;clearvars
aimpath = 'F:\EnviroVars_2020_to_2023\HYCOM_GLBy0.08_daily';
region = [262 283 18 32]; % the Gulf of Mexico
timeTick = datetime(2020,1,1):datetime(2023,12,31);
tempfolder = 'E:\Temp\HYCOM_GLBy0.08';

varList = {'temp','salt','u','v'};
URL = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0?';

for d = 1:length(timeTick)
    hourTick = timeTick(d):hours(3):timeTick(d)+day(1);
    hourTick = hourTick(1:end-1);
    nHrs = length(hourTick);
    Hourtemp = cell(nHrs,1);
    Hoursalt = cell(nHrs,1);
    HouruVel = cell(nHrs,1);
    HourvVel = cell(nHrs,1);
    Houridx = cell(nHrs,1);
    for h = 1:nHrs
        try
        D = get_hycom_online(aimpath, region, hourTick(h), varList, URL);
        if ~exist("lat",'var')
            lat = D.lat;
            lon = D.lon;
            depth = D.depth;
            time = D.time;
        end
        Hourtemp{h} = D.temp;
        Hoursalt{h} = D.salt;
        HouruVel{h} = D.u;
        HourvVel{h} = D.v;
        Houridx{h} = datestr(hourTick(h),'yyyymmddTHHMMZ');
        catch ME
            switch ME.identifier
                case 'MATLAB:imagesci:netcdf:librarydapFailure'
                    error('%s',ME.message)
                otherwise
                    fprintf('!!!-Missing data in iteration %s, skipped.\n', datestr(hourTick(h),'yyyymmddTHHMMZ'));
            end
        end
    end

    temperature = cat(4,Hourtemp{:});
    temperature(isnan(temperature)) = 0;
    dayMeanTemperature = sum(temperature, 4)./ sum(temperature ~= 0, 4);

    salinity = cat(4,Hoursalt{:});
    salinity(isnan(salinity)) = 0;
    dayMeanSalinity = sum(salinity, 4)./ sum(salinity ~= 0, 4);

    uVel = cat(4,HouruVel{:});
    uVel(isnan(uVel)) = 0;
    dayMeanUVel = sum(uVel, 4)./ sum(uVel ~= 0, 4);

    vVel = cat(4,HourvVel{:});
    vVel(isnan(vVel)) = 0;
    dayMeanVVel = sum(vVel, 4)./ sum(vVel ~= 0, 4);

    matFile = ['hycom_glby_93_',datestr(timeTick(d),'yyyymmdd'),'_daily.mat'];
    save(fullfile(aimpath,matFile),"dayMeanSalinity","dayMeanTemperature","dayMeanUVel","dayMeanVVel","lat","lon","depth","time","Houridx")
    fprintf('Saved %s\n',matFile)
    clear lat
end


