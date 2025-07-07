%location of the instrument
% lat = -(60+52.135/60);      %negative is S
% lon = 360-(56+0.628/60);    %I think the range is 0-360;

lat = -(61+15.112/60);  % EI: -(60+53.214/60);  SSI: -(61+27.469/60);  EIE: -(61+15.112/60);
lon = 360-(53+29.006/60);  % EI: 360-(55+57.238/60);  SSI: 360-(57+56.515/60)  EIE: 360-(53+29.006/60);

startDepl = '2016-02-04';  %EI: 2014-03-05;  SSI: 2015-02-10;  EIE: 2016-02-04
endDepl = '2016-12-02';  %EI: 2014-07-17;  SSI: 2015-01-29; EIE: 2016-12-02

% pull seaice data
% dataice = hdfread('K:\JST papers\Conferences\2020\Ocean Sciences (San Diego)\sea ice\asi-AMSR2-s3125-20160901-v5.4.hdf','ASI Ice Concentration');
datalat = hdfread('\\frosty.ucsd.edu\GOOGLE_DRIVE_BU\Antarctic_HARPs\Sea ice\sea ice\asi-AMSR2-s3125-20140305-v5.4.hdf','Latitudes');

datalat = hdfread('J:\antarctica\sea ice\LongitudeLatitudeGrid-s3125-Antarctic3125.hdf','Latitudes');
datalong = hdfread('J:\antarctica\sea ice\LongitudeLatitudeGrid-s3125-Antarctic3125.hdf','Longitudes');

detfn = 'asi.*.hdf';
iceDir = 'J:\antarctica\sea ice';

%min difference calculation
posdiff = sqrt((datalat-lat).^2+(datalong-lon).^2);
[x,i] = min(posdiff);
[y,j] = min(x);

%closest location is this one...
siteLat = datalat(i(j),j);
siteLong = datalong(i(j),j);

fileList = cellstr(ls(iceDir));
matchIdx = find(~cellfun(@isempty,regexp(fileList,detfn))>0);

iceData= nan(length(matchIdx),2);

for f = 1:length(matchIdx)
    % read data
    file = fileList{matchIdx(f)};
    dataiceDir = fullfile(iceDir,file);
    dataice = hdfread(dataiceDir,'ASI Ice Concentration');
    
    % get date from file
%     datetext = file(end-16:end-9);  % EI
    datetext = file(end-14:end-7);  % SSI, EIE
    
    dateNum = datenum(datetext,'yyyymmdd');
    
    iceData(f,1) = dateNum; %m2xdate(dateNum);
    iceData(f,2) = dataice(i(j),j);
end

% keep data within HARP effort range
startNum = datenum(startDepl,'yyyy-mm-dd');
endNum = datenum(endDepl,'yyyy-mm-dd');
effortTimes = find(startNum<=iceData(:,1)& iceData(:,1)<endNum);
iceData = iceData(effortTimes,:); 

figure,plot(iceData(:,1),iceData(:,2)); 

% in datetime format (easier to plot)
% dt = datetime(iceData(:,1),'ConvertFrom','datenum');
% figure,plot(dt,iceData(:,2)); 