clearvars
%close all

% Folder where to find data and store table
getFullPath = @(folder) fullfile('D:', folder);
envfolder = getFullPath('\EnviroVars_2020_to_2024\HYCOM_GLBy0.08_daily');
outfolder = getFullPath('\Results\BRT_2020_2024\Timeseries_Plots');
sitesFile = getFullPath('\Nominal_locations_sites.csv');

% Specify actions and constants:
% 1) interpolate data between depths, more accurate depth with 26sigma
stepSize = 1; % regular step size to interpolate depth (1meter)
% 2) if you wish to plot the data per file
plotting = false; % true - plot, false - no plotting
% 3) define radius in km around site to average data
radius = 8;
% 4) define specific string to match files for reading
dataFileMatchString = '\*_daily.mat';
% 5) define missing value to match files for reading
missValue = 1000; %based on netcdf files -30000, but averaged changed

% Read mat files
list = dir(fullfile(envfolder,dataFileMatchString));
list = list(~endsWith({list.name},'_surface_daily.mat'));
fileList = {list.name}';
Nfiles = length(fileList);

% Read file with site coordinates
sites = readtable(sitesFile);
sites = convertvars(sites,"SiteName",'categorical');

% Define km radius for a square around the site
deg = km2deg(radius,'earth'); %convert km radius to deg

% Read each file in folder to extract temperature and salinity.
% Crop data around the site to compute density of water, and limit to
% 26sigma to get the depth where it occurs.
Nfiles =1613;
temperature_diagram = nan(Nfiles,40);
salinity_diagram = nan(Nfiles,40);
dens_diagram = nan(Nfiles,40);
KE_diagram = nan(Nfiles,40);
vorticity_diagram = nan(Nfiles,40);
divergence_diagram = nan(Nfiles,40);
allTimes = NaT(Nfiles,1);

r=1;
fixedDim = 'Area';
TimeAxis = 1;
s = 8;%:size(sites,1);
for f = 1:Nfiles%245:(245+Nfiles)
    filename = (fullfile(envfolder,fileList{f}));
    fprintf('Reading file %s\n',filename)

    % Read latitude, longitude, depths and time from current file
    data = load(filename);
    lat = data.lat;
    lon = data.lon;
    depths = data.depth;

    % populate matrices
    allTimes(r) = data.time;

    % Read temperature, salinity, u, and v from current file
    temperature = data.dayMeanTemperature;
    salinity = data.dayMeanSalinity;
    uVel = data.dayMeanUVel;
    vVel = data.dayMeanVVel;

    %%% Check for values that are missing and change them to nan
    temperature(abs(temperature) >= missValue) = nan;
    salinity(abs(salinity) >= missValue) = nan;
    uVel(abs(uVel) >= missValue) = nan;
    vVel(abs(vVel) >= missValue) = nan;

    if sum(abs(temperature) >= missValue,"All")
        fprintf('Temperature data with %d values outside valid range\n',sum(abs(temperature) >= missValue,"All"))
    end
    if sum(abs(salinity) >= missValue,"All")
        fprintf('Salinity data with %d values outside valid range\n',sum(abs(salinity) >= missValue,"All"))
    end
    if sum(abs(uVel) >= missValue,"All")
        fprintf('u velocity data with %d values outside valid range\n',sum(uNans,"All"))
    end
    if sum(abs(vVel) >= missValue,"All")
        fprintf('v velocity data with %d values outside valid range\n',sum(vNans,"All"))
    end

    %%% Convert to Absolute Salinity and Conservative Temperature

    % build longitude/latitude grids (263×351)
    [LonP,LatP,Height] = ndgrid(lon, lat, -depths);

    % convert depth to pressure (dbar) to estimate water density
    p = gsw_p_from_z(Height, LatP);

    % calculate Absolute Salinity (SA), and Conservative Temperature (CT)
    SA = gsw_SA_from_SP(salinity,  p, LonP, LatP); % g kg⁻¹
    CT = gsw_CT_from_t (SA,temperature, p); % °C

    % 4.  Water density
    rho = gsw_rho(SA, CT, p)-1000; % in-situ water density (ρ)  (kg m⁻³) 

    %%% Compute vorticity and divergence

    % Compute the derivatives of velocity components with respect to
    % longitude (x) and latitude (y). Make sure that uVel and vVel are in
    % lon x lat size, so X x Y size.

    % Create matrix and shift them in the x-axis for longitude and y-axis
    % for latitude, so we can compute dx and dy - the distance between two
    % in the same axis.
    lonM = repmat(lon,1,length(lat));
    latM = repmat(lat',length(lon),1);

    lonM2 = circshift(lonM,1,1);
    latM2 = circshift(latM,1,2);

    % Distance between points required for estimating relative vorticity
    dx = deg2km(distance(latM,lonM,latM,lonM2))*1000;
    dy = deg2km(distance(latM,lonM,latM2,lonM))*1000;

    % Compute kinetic energy
    KE = 1/2*(uVel.^2 + vVel.^2);

    % Estimate difference between points in the x axis for u velocities
    % and in the y axis for v velocities. Will take the central
    % difference (with gradient), but for boundary points will take the
    % difference to the next point for handling nans (specifically at
    % depth).
    % central difference for interior points in the respective axis
    [du_dx, ~] = gradient(uVel);
    [~, dv_dy] = gradient(vVel);

    % Compute forward differences for boundary handling
    du_dx_boundary = diff(uVel, 1, 1);
    dv_dy_boundary = diff(vVel, 1, 2);
    du_dx_boundary = [du_dx_boundary; nan(1, size(uVel, 2), size(uVel, 3))];
    dv_dy_boundary = [dv_dy_boundary, nan(size(vVel, 1), 1, size(vVel, 3))];

    % Replace NaNs in du_dx and dv_dy with corresponding values from
    % boundary arrays
    du_dx(isnan(du_dx)) = du_dx_boundary(isnan(du_dx));
    dv_dy(isnan(dv_dy)) = dv_dy_boundary(isnan(dv_dy));

    % compute relative vorticity: dv/dx-du/dy
    vorticity = dv_dy./dx - du_dx./dy;
    % compute horizontal divergence: du/dx + dv/dy
    divergence = du_dx./dx + dv_dy./dy;

    switch fixedDim

        case 'Near' % get the nearest data to the site
            [~, lat_idx] = min(abs(lat - sites.Latitude(s)));
            [~, lon_idx] = min(abs(lon - wrapTo360(sites.Longitude(s))));

            lat_mask = (lat == lat(lat_idx));
            lon_mask = (lon == lon(lon_idx));

            temperature_cropped = squeeze(temperature(lon_mask, lat_mask,:));
            salinity_cropped = squeeze(salinity(lon_mask, lat_mask,:));
            KE_cropped = squeeze(KE(lon_mask, lat_mask,:));
            vorticity_cropped = squeeze(vorticity(lon_mask, lat_mask,:));
            divergence_cropped = squeeze(divergence(lon_mask, lat_mask,:));
            CT_cropped = squeeze(CT(lon_mask, lat_mask,:));
            SA_cropped = squeeze(SA(lon_mask, lat_mask,:));
            rho_cropped = squeeze(rho(lon_mask, lat_mask,:));

            % plot with reference density at pressure 0
            % figure,gsw_SA_CT_plot(SA_cropped,CT_cropped,0,'\itS\rm_A - \Theta plot') 


        case 'Area' % Compute the average of the data within the radius of the site
            topLat = sites.Latitude(s) + deg;
            bottomLat = sites.Latitude(s) - deg;
            leftLon = sites.Longitude(s) - deg;
            rightLon = sites.Longitude(s) + deg;

            crop_lat_range = [topLat bottomLat];
            crop_lon_range = wrapTo360([leftLon rightLon]);

            lat_mask = (lat >= crop_lat_range(2)) & (lat <= crop_lat_range(1));
            lon_mask = (lon >= crop_lon_range(1)) & (lon <= crop_lon_range(2));

            % Apply the logical masks to obtain the cropped data
            temperature_area = temperature(lon_mask, lat_mask,:);
            salinity_area = salinity(lon_mask, lat_mask,:);
            rho_area = rho(lon_mask, lat_mask,:);
            KE_area = KE(lon_mask, lat_mask,:);
            vorticity_area = vorticity(lon_mask, lat_mask,:);
            divergence_area = divergence(lon_mask, lat_mask,:);
            SA_area = SA(lon_mask, lat_mask,:);
            CT_area = CT(lon_mask, lat_mask,:);

            % Average data within site limits to obtain one value per depth
            avg_temperature = mean(temperature_area,[1 2],"omitnan");
            temperature_cropped = reshape(avg_temperature,[size(avg_temperature,3),1]); % 3D matrix 1x1x40 to vector 40x1

            avg_salinity = mean(salinity_area,[1 2],"omitnan");
            salinity_cropped = reshape(avg_salinity,[size(avg_salinity,3),1]); % 3D matrix 1x1x40 to vector 40x1

            avg_rho = mean(rho_area,[1 2],"omitnan");
            rho_cropped = reshape(avg_rho,[size(avg_rho,3),1]); % 3D matrix 1x1x40 to vector 40x1

            avg_KE = mean(KE_area,[1 2],"omitnan");
            KE_cropped = reshape(avg_KE,[size(avg_KE,3),1]); % 3D matrix 1x1x40 to vector 40x1

            avg_vorticity = mean(vorticity_area,[1 2],"omitnan");
            vorticity_cropped = reshape(avg_vorticity,[size(avg_vorticity,3),1]); % 3D matrix 1x1x40 to vector 40x1

            avg_divergence = mean(divergence_area,[1 2],"omitnan");
            divergence_cropped = reshape(avg_divergence,[size(avg_divergence,3),1]); % 3D matrix 1x1x40 to vector 40x1

            avg_CT = mean(CT_area,[1 2],"omitnan");
            CT_cropped = reshape(avg_CT,[size(avg_CT,3),1]); % 3D matrix 1x1x40 to vector 40x1

            avg_SA = mean(SA_area,[1 2],"omitnan");
            SA_cropped = reshape(avg_SA,[size(avg_SA,3),1]); % 3D matrix 1x1x40 to vector 40x1
    end

        temperature_diagram(r,1:length(temperature_cropped)) = temperature_cropped(:)';
        salinity_diagram(r,1:length(salinity_cropped)) = salinity_cropped(:)';
        dens_diagram(r,1:length(rho_cropped)) = rho_cropped(:)';
        
        KE_diagram(r,1:length(KE_cropped)) = KE_cropped(:)';
        vorticity_diagram(r,1:length(vorticity_cropped)) = vorticity_cropped(:)';
        divergence_diagram(r,1:length(divergence_cropped)) = divergence_cropped(:)';

        CT_diagram(r,1:length(CT_cropped)) = CT_cropped(:)';
        SA_diagram(r,1:length(SA_cropped)) = SA_cropped(:)';


    r=r+1;
end


dayEnv = timetable(allTimes,temperature_diagram,salinity_diagram,dens_diagram,KE_diagram,vorticity_diagram,divergence_diagram,SA_diagram,CT_diagram);
dayEnv.Properties.DimensionNames{'allTimes'} = 'Time';

% Get species densities to combine to one table
tablesPath = getFullPath('\Results\populated_tables');
load(fullfile(tablesPath,'All_sites_2020_2024_densities_and_CV.mat'))
dayDens = dayDens(dayDens.Site == sites.SiteName(s), :); % filter to specific site

combined = synchronize(dayEnv, dayDens, 'intersection');
combined.dens_diagram = combined.dens_diagram;
combined.year = year(combined.Time);


outmatfile = sprintf('Site_%s_2020_2024_SpDens_HYCOM_avg_%dkm_radius_CT_SA.mat',sites.SiteName(s),radius);
save(fullfile(outfolder,outmatfile),'combined','depths')



%%%% Plotting

%%%%% 5-day smoothing
% Temperature, Salinity, Density, KE, Vorticity

% Specific color palette for kinetic energy
colorKE = [1	1	1
    0.900000000000000	0.944444444444444	0.944444444444444
    0.800000000000000	0.888888888888889	0.888888888888889
    0.700000000000000	0.833333333333333	0.833333333333333
    0.600000000000000	0.777777777777778	0.777777777777778
    0.500000000000000	0.722222222222222	0.722222222222222
    0.400000000000000	0.666666666666667	0.666666666666667
    0.300000000000000	0.611111111111111	0.611111111111111
    0.200000000000000	0.555555555555556	0.555555555555556
    0.100000000000000	0.500000000000000	0.500000000000000
    0.107502664818649	0.490164566371003	0.465646767823488
    0.165903166962839	0.513789544168112	0.458655372033719
    0.224303669107030	0.537414521965222	0.451663976243951
    0.282704171251220	0.561039499762332	0.444672580454182
    0.341104673395410	0.584664477559442	0.437681184664414
    0.399505175539600	0.608289455356551	0.430689788874646
    0.457905677683790	0.631914433153661	0.423698393084877
    0.516306179827980	0.655539410950771	0.416706997295109
    0.574706681972171	0.679164388747881	0.409715601505340
    0.635441729573654	0.683957813630852	0.368731275428919
    0.696176777175137	0.688751238513822	0.327746949352498
    0.722402524220031	0.685914192712434	0.328696514477442
    0.748628271264926	0.683077146911046	0.329646079602385
    0.774854018309821	0.680240101109658	0.330595644727329
    0.801079765354715	0.677403055308270	0.331545209852273
    0.827305512399610	0.674566009506882	0.332494774977216
    0.853531259444505	0.671728963705495	0.333444340102160
    0.859335756028145	0.664981590393718	0.375794821006860
    0.865140252611785	0.658234217081942	0.418145301911561
    0.870944749195425	0.651486843770166	0.460495782816261
    0.876749245779065	0.644739470458390	0.502846263720962
    0.882553742362705	0.637992097146614	0.545196744625662
    0.869066133348282	0.586472240198754	0.496176740551260
    0.854076486576040	0.534952196157316	0.455228518937465
    0.836399335601840	0.483998874468285	0.423549487463078
    0.815079296885781	0.434340532127166	0.401255316597207
    0.789637141880046	0.386642687691636	0.387277739793126
    0.760160645444544	0.341286359687873	0.379711630258191
    0.727025555782371	0.298393902856668	0.376595129082709
    0.690686821415406	0.257930683487536	0.376219604697688
    0.651502778140272	0.219872649535767	0.377199995318858
    0.609670032026622	0.184371393079528	0.378317256256392
    0.565221845822222	0.151950786117574	0.378284503952529
    0.518037893679901	0.123803312924417	0.375433299742202
    0.468041751020257	0.101737226058069	0.367511084450080
    0.415530934232078	0.0870820265879038	0.351892976555934
    0.361556264060448	0.0783842163701545	0.326932726816117
    0.307591897201711	0.0716212985362024	0.293220368404476
    0.254744188185617	0.0631875246166560	0.252979256517997
    0.203400246337400	0.0512571528517886	0.208537206377127];

windowSize = 5;
% compute 5-day moving averages, ignoring NaNs
combined.smMd = movmean(combined.densMd, windowSize, 'omitnan');
combined.smMe = movmean(combined.densMe, windowSize, 'omitnan');
combined.smZc = movmean(combined.densZc, windowSize, 'omitnan');
combined.smPm = movmean(combined.densPm, windowSize, 'omitnan');

combined.smT = movmean(combined.temperature_diagram, windowSize, 'omitnan');
combined.smS = movmean(combined.salinity_diagram, windowSize, 'omitnan');
combined.smD = movmean(combined.dens_diagram, windowSize, 'omitnan');
combined.smKE = movmean(combined.KE_diagram, windowSize, 'omitnan');
combined.smVort = movmean(combined.vorticity_diagram, windowSize, 'omitnan');
combined.smDiver = movmean(combined.divergence_diagram, windowSize, 'omitnan');


figure;
Nplots = 10;
tiledlayout(Nplots, 1, 'TileSpacing', 'none', 'Padding', 'none');

ax1 = nexttile;
% Plot Md density
plot(datenum(combined.Time), combined.smMd, 'k-', 'LineWidth', 1.5)
set(ax1,'XTickLabel', [])  % hide x-axis labels
ylabel('Md Density')
title(sprintf('Site %s (%d-day moving average)',SiteName,windowSize))

ax2 = nexttile;
% Plot Me density
plot(datenum(combined.Time), combined.smMe, 'k-', 'LineWidth', 1.5)
set(ax2,'XTickLabel', [])  % hide x-axis labels
ylabel('Me Density')

ax3 = nexttile;
% Plot Zc density
plot(datenum(combined.Time), combined.smZc, 'k-', 'LineWidth', 1.5)
set(ax3,'XTickLabel', [])  % hide x-axis labels
ylabel('Zc Density')

ax4 = nexttile;
% Plot Pm density
plot(datenum(combined.Time), combined.smPm, 'k-', 'LineWidth', 1.5)
set(ax4,'XTickLabel', [])  % hide x-axis labels
ylabel('Pm Density')

% Water column plots
nanCols = all(isnan(combined.dens_diagram), 1);
firstAllNaNCol = find(nanCols, 1, 'first');
Z = 1:firstAllNaNCol-1;
%Z = 27:firstAllNaNCol-1;

% Original depths
original_depths = depths(Z);

% Create pseudo-log spacing
stretched_depths = log10(original_depths + 10); % +10 to avoid log(0)

% Normalize to match the original scale for plotting
stretched_depths = stretched_depths - min(stretched_depths);
stretched_depths = stretched_depths / max(stretched_depths) * max(original_depths);

% Preallocate and store axis handles
axList = gobjects(Nplots, 1);
axList(1:4) = [ax1, ax2, ax3, ax4];

if Z(1) == 1
    yticks_to_show = [0 10 20 50 100 200 500 1000 2000 4000];
else
    yticks_to_show = [0 10 20 50 100 200 400 600 800 1000 1500 2000 3000 4000];
end
ytick_positions = interp1(original_depths, stretched_depths, yticks_to_show, 'linear', 'extrap');

for wM = 1:6
    ax = nexttile;
    axList(4 + wM) = ax;
    gLim=[];
    switch wM
        case 1
            paramProfile = combined.smD(:,Z);
            paramName = 'Water Density';
        case 2
            paramProfile = combined.smT(:,Z);
            paramName = 'Temperature';
        case 3
            paramProfile = combined.smS(:,Z);
            paramName = 'Salinity';
        case 4
            paramProfile = combined.smKE(:,Z);
            paramName = 'Kinetic Energy';
        case 5
            paramProfile = combined.smVort(:,Z);
            paramName = 'Vorticity';
        case 6
            paramProfile = combined.smDiver(:,Z);
            paramName = 'Divergence';
    end

    pcolor(datenum(combined.Time), stretched_depths, paramProfile')
    shading interp
    switch wM
        case 1
            %colormap(ax, cmocean('delta',32,'pivot',26))
            colormap(ax, flipud(colorKE))
            %gLim = [20, 28];
        case 2
            colormap(ax, cmocean('thermal',32))
            %gLim = [0, 33];
        case 3
            colormap(ax, cmocean('haline',26))
            %gLim = [32.5, 37.5];
        case 4
            colormap(ax, colorKE)
            gLim = [0 max(paramProfile,[],'all')];
        case 5
            colormap(ax, cmocean('-balance',30))
            gLim = [-max(abs(paramProfile),[],'all'),max(abs(paramProfile),[],'all')];
        case 6
            colormap(ax, cmocean('-balance',30))
            gLim = [-max(abs(paramProfile),[],'all'),max(abs(paramProfile),[],'all')];
        
    end
    if isempty(gLim)
        set(ax,'YTick', ytick_positions, 'YTickLabel', string(yticks_to_show), ...
            'YDir', 'reverse')
    else
        set(ax, 'CLim',gLim,'YTick', ytick_positions, 'YTickLabel', string(yticks_to_show), ...
            'YDir', 'reverse')
    end
    if wM==1
        hold(ax,'on')
        [C, h] = contour(datenum(combined.Time), stretched_depths, paramProfile', ...
            [22 26,27.44,27.80], 'LineColor', 'k', 'LineWidth', 1);
        clabel(C, h, 'FontSize', 8, 'Color', 'k');
        hold(ax,'off')
    end

    ylabel('Depth')
    grid(ax, 'on')
    set(ax, 'Layer', 'top')

    % Add colorbar with label
    cb = colorbar(ax);
    ylabel(cb, paramName)

    % Hide x-axis labels for all but last tile
    if wM < 6
        set(ax, 'XTickLabel', [])
    end
end

datetick(ax, 'x', 'mmm''yy', 'keeplimits', 'keepticks')
% Link all x-axes
linkaxes(axList, 'x');
linkaxes(axList(5:end), 'y');

% Identify time periods with NaNs in Zc
dNumEffort = datenum(combined.Time);
isNanEffort = isnan(combined.densZc);
nanBlocks = bwlabel(isNanEffort);
blockIDs = unique(nanBlocks(nanBlocks > 0));

% Overlay patches on first 4 axes
for a = 1:4
    ax = axList(a);
    hold(ax, 'on')
    yLimits = ylim(ax);

    for i = 1:length(blockIDs)
        idx = find(nanBlocks == blockIDs(i));
        x1 = dNumEffort(idx(1));
        x2 = dNumEffort(idx(end));
        patch(ax, [x1 x2 x2 x1], ...
            [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
            [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    end

    hold(ax, 'off')
end

set(axList, 'XLim', [min(dNumEffort), max(dNumEffort)]);

% Vector with the first day of every month that the series spans
minorDt = dateshift(combined.Time(1),'start','month'):calmonths(1):dateshift(combined.Time(end),'start','month');
majorDt = dateshift(combined.Time(1),'start','year'):calyears(1):dateshift(combined.Time(end),'start','year');

minorNums = datenum(minorDt);
majorNums = datenum(majorDt);

% Put the same monthly ticks on every axis
set(axList , 'XTick' , majorNums,'XMinorTick','on')

for ax = axList(:)'
    ax.XAxis.MinorTickValues = minorNums;
end

% Put month labels on (only need to label one axis once they’re linked)
datetick(axList(end) , 'x' , 'mmm''yy' , 'keeplimits' , 'keepticks')

