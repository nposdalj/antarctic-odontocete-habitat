%% make a map of SOCAL
% from LMB last update 10/09/24
% CMS updated for SOCAL Baleen Whales
% load in CASE grid data:
load('C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Site Map/GMRT_bathymetry/grid.mat') % XC - lon, YC - lat, RC - depth
% get lon in +/-180 vs 0-360
XC =  -(360-XC);

% bathymetry data from gmrt, extract as a grd file
% [X, Y, Z] = grdread2('D:\Other\Bathymetry\CADEMO_zoom.grd'); % load the data
% contour(X,Y,Z,'black','showtext','on'); % check that it looks right
% save('D:\Other\Bathymetry\CHNMS_NO_bathy.mat','X','Y','Z') % save if you want
% load('F:\Tracking\bathymetry\socal_new.mat');

% % read in the polygons
% s = shaperead('F:\bathymetry\cinms_py2\cinms_py.shp','UseGeoCoords',true);
%{
% load in our site lat/lon
load('G:\My Drive\Maps\socal2.mat');
v = [0,0];
% CHNMS_NO = [34.5743 -120.7133];
% EI = [(60+53.214/60)*-1 (55+57.238/60)*-1 762];
% EIE = [(61+15.112/60)*-1 (53+29.006/60)*-1 1033];
% SSI = [(61+27.469/60)*-1 (57+56.515/60)*-1 768];
% % l = 0:10:max(Z,[],'all');
H = [32.86117  -119.13516 -1282.5282];
W = [33.53973  -120.25815 -1377.8741];
N = [32.36975  -118.56458 -1298.3579];
E = [32.65345  -119.48455 -1328.9836];

% read kmz file!
% T = kmz2struct('D:\Other\CADEMO.kmz');
% read in polygons
% s = shaperead('D:\Other\CA_cst3nm.shp','UseGeoCoords',true);


% define the colormap
% a = cmocean('ice')
a = m_colmap('blues'); % for colormap: https://urldefense.com/v3/__https://github.com/g2e/m_map/blob/master/m_colmap.m__;!!Mih3wA!CmpmXt6-4iAp26KsnI7jGzbJ2UPd4vzKZ66z1YZcF3Q-6bwxmDFfaRYu3AIGRvXd73BiLKSFc4472kt2SA$ 
% a = a(100:200,:)
a = vertcat(a,[0.8020    0.8020    0.8020]); % add grey for land
% a = a(25:end,:);

figure
contourf(X,Y,Z,50,'edgecolor','none') %10) %,'edgecolor','black') %'showtext','off'
colormap(a)
caxis([-2000 0])
% hold on
% geoshow(s,'facecolor','none','edgecolor',[1.0000    0.4118    0.1608],'LineWidth',3)
hold on
contourf(X,Y,Z,v,'k')
% c = colorbar
% ylabel(c,'Depth (m)')
% xticklabels({'121.5°W','121°W','120.5°W','120°W','119.5°W','119°W'})
% yticklabels({'33.5°N', '34°N','34.5°N','35°N','35.5°N'})
% scatter(-119.6982,34.4208,20,'o','filled','black')
% scatter(-120.8500,35.3659,20,'o','filled','black')
% scatter(-120.6599,35.2825,20,'o','filled','black')
% scatter(-120.4716,34.4486,20,'o','filled','black')

xlim([min(X) max(X)])
ylim([min(Y) max(Y)])

scatter(XC,YC,7,'o','filled','yellow')
scatter(W(2),W(1),80,'o','filled','red')
scatter(H(2),H(1),80,'o','filled','red')
scatter(E(2),E(1),80,'o','filled','red')
scatter(N(2),N(1),80,'o','filled','red')
yticks([32 33 34 35])
yticklabels({'32°N', '33°N','34°N','35°N'})
ytickangle(90)
xticks([-120, -119, -118])
xticklabels({'120°W','119°W','118°W'})
% legend(b,{'SOCAL H'})
% add a point for landmarks
scatter(-118.2437,34.0722,20,'o','filled','black')
scatter(-119.6982,34.5208,20,'o','filled','black')
text(-118.2437,34.0722,'Los Angeles','HorizontalAlignment','right','fontsize',12,'fontweight','bold')
text(W(2),W(1)+0.07,'W','HorizontalAlignment','center','fontsize',12,'fontweight','bold')
text(H(2),H(1)+0.07,'H','HorizontalAlignment','center','fontsize',12,'fontweight','bold')
text(E(2),E(1)+0.07,'E','HorizontalAlignment','center','fontsize',12,'fontweight','bold')
text(N(2),N(1)+0.07,'N','HorizontalAlignment','center','fontsize',12,'fontweight','bold')
text(-118.4981,32.9029,'San Clemente','HorizontalAlignment','center','fontsize',8,'fontweight','bold','FontName','Georgia')
text(-119.4992,33.2465,'San Nicholas','HorizontalAlignment','center','fontsize',8,'fontweight','bold','FontName','Georgia')
text(-118.4163,33.3879,'Catalina','HorizontalAlignment','center','fontsize',8,'fontweight','bold','FontName','Georgia')
text(-119.7658,34.0232,'Santa Cruz','HorizontalAlignment','center','fontsize',8,'fontweight','bold','FontName','Georgia')
text(-120.0896,33.9773,'Santa Rosa','HorizontalAlignment','center','fontsize',8,'fontweight','bold','FontName','Georgia')
text(-120.3724,34.0376,'San Miguel','HorizontalAlignment','center','fontsize',7)
text(-119.6982,34.5208,'Santa Barbara','HorizontalAlignment','right','fontsize',12,'fontweight','bold')

%% zoomed in map with turbine locations

[X, Y, Z] = grdread2('D:\Other\Bathymetry\CADEMO_zoom.grd'); % load the data
contour(X,Y,Z,'black','showtext','on'); % check that it looks right

% read kmz file!
T = kmz2struct('D:\Other\CADEMO.kmz');

figure
contour(X,Y,Z,'showtext','on') % ,'edgecolor','none')
contourf(X,Y,Z,'showtext','on') % ,'edgecolor','none')
colormap(a)
caxis([-125 0])
% hold on
% geoshow(s,'facecolor','none','edgecolor','red')
hold on
% contourf(X,Y,Z,v,'k','ShowText','on')
% c = colorbar
% ylabel(c,'Depth (m)')
scatter(CHNMS_NO(2),CHNMS_NO(1),50,'o','filled','red')
% xticklabels({'121.5°W','121°W','120.5°W','120°W','119.5°W','119°W'})
% yticklabels({'33.5°N', '34°N','34.5°N','35°N','35.5°N'})
% scatter(-119.6982,34.4208,20,'o','filled','black')
% scatter(-120.8500,35.3659,20,'o','filled','black')
% scatter(-120.6599,35.2825,20,'o','filled','black')
% scatter(-120.4716,34.4486,20,'o','filled','black')
% add wind turbine sites
scatter(T(1).Lon,T(1).Lat,'^','filled','yellow')
scatter(T(2).Lon,T(2).Lat,'^','filled','yellow')
scatter(T(3).Lon,T(3).Lat,'^','filled','yellow')
scatter(T(4).Lon,T(4).Lat,'^','filled','yellow')

%}
%% site map for antarctic region!
% bathymetry data from gmrt, extract as a grd file
[X, Y, Z] = grdread2('C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Site Map/GMRT_bathymetry/GMRTv4_3_1_20250703topo.grd'); % load the data
%[A,R] = readgeoraster("F:\Antarctica\bathymetry\antTif.tif","outputtype","double",'CoordinateSystemType','geographic');
latlim = [-62 -60];
lonlim = [-59 -52];

[lonGrid, latGrid] = meshgrid(X, Y);  % grdread2 gives X = lon, Y = lat vectors
R = georasterref('RasterSize', size(Z), ...
                 'Latlim', [min(latGrid(:)) max(latGrid(:))], ...
                 'Lonlim', [min(lonGrid(:)) max(lonGrid(:))], ...
                 'RasterInterpretation', 'cells');

% load('F:\bathymetry\socal_new.mat');
% v = [0,0];
EI = [(60+53.214/60)*-1 (55+57.238/60)*-1 762];
CI = [(61+15.112/60)*-1 (53+29.006/60)*-1 1033];
KGI = [(61+27.469/60)*-1 (57+56.515/60)*-1 768];

%a = cmocean('deep')
a = m_colmap('blue',256);
a = vertcat(a, [0.8020 0.8020 0.8020]); % add grey

figure
worldmap(latlim,lonlim)
%geoshow(A,R,'displaytype','surface')
geoshow(Z, R, 'DisplayType', 'surface');

colormap(a)
caxis([-2000 0])
cb = colorbar;
ylabel(cb, 'Depth (m)')

% plot(EI(2),EI(1),'markersize',150,'markerfacecolor','red')
setm(gca, 'MLabelParallel', -62, 'MeridianLabel', 'on', 'MLabelLocation', 2)
%gs = scaleruler('Units','km','RulerStyle','patches','MajorTick',50,'MinorTick',5);
%set(gs, 'XLoc', -56.5, 'YLoc', -61.9);  % Lat/Lon coordinates for placement

geoshow(EI(1),EI(2), 'DisplayType', 'Point', 'Marker', 'o', 'Color', 'red','markerfacecolor','red');
geoshow(CI(1),CI(2), 'DisplayType', 'Point', 'marker','o','Color', 'red','markerfacecolor','red');
geoshow(KGI(1),KGI(2), 'DisplayType', 'Point','marker','o', 'Color', 'red','markerfacecolor','red');
legend off
