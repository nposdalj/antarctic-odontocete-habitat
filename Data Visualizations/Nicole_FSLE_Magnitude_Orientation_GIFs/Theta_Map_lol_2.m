indir = 'E:\New_FSLE_240612';

% output directory for figures
outdir = 'E:\New_FSLE_240612';

%position CCE1
latCCE1 = 32.8986;
lonCCE1 = 242.3902;

% size of plot in degrees
s = 2;
% 6 for larger plot, 2 for smaller plot

% number of pixel n (uneven number) for fsle value at site
%these pixels are the exact pixels used in FSLE grid cells to get data from
% kept seperate lat and lon n to change suqrae size for each range
lon_n = 11;
lat_n =9;

fsle_mean = ones(1,1)*NaN;
theta_mean = ones(1,1)*NaN;
fsle_std = ones(1,1)*NaN;
theta_std = ones(1,1)*NaN;
matTime = ones(1,1)*NaN;




%read in data from nc file
    fsleMaxData  = ncread('E:\New_FSLE_240612\dt_global_allsat_madt_fsle_20170216_20180704.nc','fsle_max');
    longitude = ncread('E:\New_FSLE_240612\dt_global_allsat_madt_fsle_20170216_20180704.nc','lon');
    latitude = ncread('E:\New_FSLE_240612\dt_global_allsat_madt_fsle_20170216_20180704.nc','lat');
    thetaMaxData = ncread('E:\New_FSLE_240612\dt_global_allsat_madt_fsle_20170216_20180704.nc', 'theta_max');

    %date
    time = ncread('E:\New_FSLE_240612\dt_global_allsat_madt_fsle_20170216_20180704.nc','time');
    time = double(time);
    adj = datenum([1950 1 1 0 0 0]);
    matTime = time + adj;
    
    % n x n pixel square around mooring location
    [~, lonCCE1Idx] = min(abs(longitude-lonCCE1));
    [~, latCCE1Idx] = min(abs(latitude-latCCE1));
    
    lonCCE1min = longitude(lonCCE1Idx-(lon_n-1)/2);
    lonCCE1max = longitude(lonCCE1Idx+(lon_n-1)/2);
    latCCE1min = latitude(latCCE1Idx-(lat_n-1)/2);
    latCCE1max = latitude(latCCE1Idx+(lat_n-1)/2);
    
 
    fsle_CCE1 = [];
    fsle_CCE1 = fsleMaxData(lonCCE1Idx-(lon_n-1)/2:lonCCE1Idx+(lon_n-1)/2, ...
        latCCE1Idx-(lat_n-1)/2:latCCE1Idx+(lat_n-1)/2);
% calculate mean/std/absmean over theta near site
    fsle_mean = mean(fsle_CCE1(:));
    fsle_std = std(fsle_CCE1(:));

    
    theta_CCE1 = [];
    theta_CCE1 = thetaMaxData(lonCCE1Idx-(lon_n-1)/2:lonCCE1Idx+(lon_n-1)/2, ...
        latCCE1Idx-(lat_n-1)/2:latCCE1Idx+(lat_n-1)/2);
    % calculate mean/std over theta near site
    theta_mean = mean(theta_CCE1(:));
    theta_std = std(theta_CCE1(:));

    
    % matrix to plot the square onto map
    lonSqu = [lonCCE1min lonCCE1max lonCCE1max lonCCE1min lonCCE1min];
    latSqu = [latCCE1max latCCE1max latCCE1min latCCE1min latCCE1max];

    % s deg square around mooring site for plotting
    long_min = lonCCE1-s/2;
    long_max = lonCCE1+s/2;
    lat_min = latCCE1-s/2;
    lat_max = latCCE1+s/2;

    x = [];
    lon = [];
    lat = [];
    [~, lonMinIdx] = min(abs(longitude-long_min));
    [~, lonMaxIdx] = min(abs(longitude-long_max));
    [~, latMinIdx] = min(abs(latitude-lat_min));
    [~, latMaxIdx] = min(abs(latitude-lat_max));
    
    % Find FSLE value for plotting
    fsle_Value = [];
    fsle_Value = fsleMaxData(lonMinIdx:lonMaxIdx, latMinIdx:latMaxIdx);
    
    % Find theta value for plotting
    theta_Value = [];
    theta_Value = thetaMaxData(lonMinIdx:lonMaxIdx, latMinIdx:latMaxIdx);
    
    
    % load bathymetry data
   % bath_lat = ncread("P:\FSLE\GEBCO_14_Feb_2022_393eed1ebda3\gebco_2021_n33.8986_s31.8986_w-118.6098_e-116.42.nc", "lat");
   % bath_lon = ncread("P:\FSLE\GEBCO_14_Feb_2022_393eed1ebda3\gebco_2021_n33.8986_s31.8986_w-118.6098_e-116.42.nc", "lon");
   % bath_elev = ncread("P:\FSLE\GEBCO_14_Feb_2022_393eed1ebda3\gebco_2021_n33.8986_s31.8986_w-118.6098_e-116.42.nc", "elevation");
    %bath_long360 = bath_lon + 360;
    %contourf(bath_long360, bath_lat, bath_elev')
   % colormap('white');
 
    
      

    %plot theta map
    f2 = subplot(1, 2, 2)
    x = [long_min:.04:long_max];
    y = [lat_min:.04:lat_max];
    R = rot90(fliplr(theta_Value));
    h1 = pcolor(x,y,R)
    crameri('corko')
    shading flat
    %cmocean('phase', 5)
    % colormap('phase', 5);
   % c = colorbar;
   % caxis([-89.9 89.9])
   %colormap('dir360colormap');
   caxis([-90 90])
   c = colorbar('Ticks',[-90,-45,0,45,90],...
      'TickLabels',{'-90 E <-> W', '-45 SW <-> NE', '0 N <-> S', '45 NW <-> SE', '90 E <-> W'});
    % c = colorbar('XTickLabel', {'', '-45  NW<->SE', '0  E<->W', '45  SW<->NE', '' }, 'XTick', -90:45:90);
   % c = colorbar('XTickLabel', {'-90  N<->S', '-45  NW<->SE', '0  E<->W', '45  SW<->NE', '90  N<->S'}, 'XTick', -90:45:90);
    %cmocean('-deep') %Comment out if you don't have cmocean
    % colormap(flipud(parula))
    ylabel(c,'FSLV Orientation (°)/FSLE Filament Orientation [per day]')
    %load coast
    hold (f2, 'on')
    plot(lonSqu,latSqu,'k','LineWidth',2)
    title(datestr(matTime))
    xlabel('Longitude')
    ylabel('Latitude')
    hold (f2, 'off')
    set(gcf, 'Position',  [100, 100, 1500, 500])
    
    %save figure to output directory
    file = strrep('dt_global_allsat_madt_fsle_20170216_20180704_corko.nc','nc','png');
    outFile = fullfile(outdir,file);
    saveas(gcf,outFile,'png')
    
    
    
    
    
    
    
    
    
    
    
    
%     close all
%     clear all
%     
%     indir = 'P:\New_FSLE_240612';
% 
% % output directory for figures
% outdir = 'P:\New_FSLE_240612';
% 
% %position CCE1
% latCCE1 = 32.8986;
% lonCCE1 = 242.3902;
% 
% % size of plot in degrees
% s = 2;
% % 6 for larger plot, 2 for smaller plot
% 
% % number of pixel n (uneven number) for fsle value at site
% %these pixels are the exact pixels used in FSLE grid cells to get data from
% % kept seperate lat and lon n to change suqrae size for each range
% lon_n = 11;
% lat_n =9;
% 
% fsle_mean = ones(1,1)*NaN;
% theta_mean = ones(1,1)*NaN;
% fsle_std = ones(1,1)*NaN;
% theta_std = ones(1,1)*NaN;
% matTime = ones(1,1)*NaN;
% 
% 
% 
% 
% %read in data from nc file
%     fsleMaxData  = ncread('P:\New_FSLE_240612\dt_global_allsat_madt_fsle_20170418_20180704.nc','fsle_max');
%     longitude = ncread('P:\New_FSLE_240612\dt_global_allsat_madt_fsle_20170418_20180704.nc','lon');
%     latitude = ncread('P:\New_FSLE_240612\dt_global_allsat_madt_fsle_20170418_20180704.nc','lat');
%     thetaMaxData = ncread('P:\New_FSLE_240612\dt_global_allsat_madt_fsle_20170418_20180704.nc', 'theta_max');
% 
%     %date
%     time = ncread('P:\New_FSLE_240612\dt_global_allsat_madt_fsle_20170418_20180704.nc','time');
%     time = double(time);
%     adj = datenum([1950 1 1 0 0 0]);
%     matTime = time + adj;
%     
%     % n x n pixel square around mooring location
%     [~, lonCCE1Idx] = min(abs(longitude-lonCCE1));
%     [~, latCCE1Idx] = min(abs(latitude-latCCE1));
%     
%     lonCCE1min = longitude(lonCCE1Idx-(lon_n-1)/2);
%     lonCCE1max = longitude(lonCCE1Idx+(lon_n-1)/2);
%     latCCE1min = latitude(latCCE1Idx-(lat_n-1)/2);
%     latCCE1max = latitude(latCCE1Idx+(lat_n-1)/2);
%     
%  
%     fsle_CCE1 = [];
%     fsle_CCE1 = fsleMaxData(lonCCE1Idx-(lon_n-1)/2:lonCCE1Idx+(lon_n-1)/2, ...
%         latCCE1Idx-(lat_n-1)/2:latCCE1Idx+(lat_n-1)/2);
% % calculate mean/std/absmean over theta near site
%     fsle_mean = mean(fsle_CCE1(:));
%     fsle_std = std(fsle_CCE1(:));
% 
%     
%     theta_CCE1 = [];
%     theta_CCE1 = thetaMaxData(lonCCE1Idx-(lon_n-1)/2:lonCCE1Idx+(lon_n-1)/2, ...
%         latCCE1Idx-(lat_n-1)/2:latCCE1Idx+(lat_n-1)/2);
%     % calculate mean/std over theta near site
%     theta_mean = mean(theta_CCE1(:));
%     theta_std = std(theta_CCE1(:));
% 
%     
%     % matrix to plot the square onto map
%     lonSqu = [lonCCE1min lonCCE1max lonCCE1max lonCCE1min lonCCE1min];
%     latSqu = [latCCE1max latCCE1max latCCE1min latCCE1min latCCE1max];
% 
%     % s deg square around mooring site for plotting
%     long_min = lonCCE1-s/2;
%     long_max = lonCCE1+s/2;
%     lat_min = latCCE1-s/2;
%     lat_max = latCCE1+s/2;
% 
%     x = [];
%     lon = [];
%     lat = [];
%     [~, lonMinIdx] = min(abs(longitude-long_min));
%     [~, lonMaxIdx] = min(abs(longitude-long_max));
%     [~, latMinIdx] = min(abs(latitude-lat_min));
%     [~, latMaxIdx] = min(abs(latitude-lat_max));
%     
%     % Find FSLE value for plotting
%     fsle_Value = [];
%     fsle_Value = fsleMaxData(lonMinIdx:lonMaxIdx, latMinIdx:latMaxIdx);
%     
%     % Find theta value for plotting
%     theta_Value = [];
%     theta_Value = thetaMaxData(lonMinIdx:lonMaxIdx, latMinIdx:latMaxIdx);
%     
%     
%     % load bathymetry data
%    % bath_lat = ncread("P:\FSLE\GEBCO_14_Feb_2022_393eed1ebda3\gebco_2021_n33.8986_s31.8986_w-118.6098_e-116.42.nc", "lat");
%    % bath_lon = ncread("P:\FSLE\GEBCO_14_Feb_2022_393eed1ebda3\gebco_2021_n33.8986_s31.8986_w-118.6098_e-116.42.nc", "lon");
%    % bath_elev = ncread("P:\FSLE\GEBCO_14_Feb_2022_393eed1ebda3\gebco_2021_n33.8986_s31.8986_w-118.6098_e-116.42.nc", "elevation");
%     %bath_long360 = bath_lon + 360;
%     %contourf(bath_long360, bath_lat, bath_elev')
%    % colormap('white');
%  
%     
%       
% 
%     %plot theta map
%     f2 = subplot(1, 2, 2)
%     x = [long_min:.04:long_max];
%     y = [lat_min:.04:lat_max];
%     R = rot90(fliplr(theta_Value));
%     h1 = pcolor(x,y,R)
%     crameri('viko')
%     shading flat
%     %cmocean('phase', 5)
%     % colormap('phase', 5);
%    % c = colorbar;
%    % caxis([-89.9 89.9])
%    %colormap('dir360colormap');
%    caxis([-90 90])
%    c = colorbar('Ticks',[-90,-45,0,45,90],...
%       'TickLabels',{'-90 E <-> W', '-45 SW <-> NE', '0 N <-> S', '45 NW <-> SE', '90 E <-> W'});
%     % c = colorbar('XTickLabel', {'', '-45  NW<->SE', '0  E<->W', '45  SW<->NE', '' }, 'XTick', -90:45:90);
%    % c = colorbar('XTickLabel', {'-90  N<->S', '-45  NW<->SE', '0  E<->W', '45  SW<->NE', '90  N<->S'}, 'XTick', -90:45:90);
%     %cmocean('-deep') %Comment out if you don't have cmocean
%     % colormap(flipud(parula))
%     ylabel(c,'FSLV Orientation (°)/FSLE Filament Orientation [per day]')
%     %load coast
%     hold (f2, 'on')
%     plot(lonSqu,latSqu,'k','LineWidth',2)
%     title(datestr(matTime))
%     xlabel('Longitude')
%     ylabel('Latitude')
%     hold (f2, 'off')
%     set(gcf, 'Position',  [100, 100, 1500, 500])
%     
%     %save figure to output directory
%     file = strrep('dt_global_allsat_madt_fsle_20170418_20180704.nc','nc','png');
%     outFile = fullfile(outdir,file);
%     saveas(gcf,outFile,'png') 