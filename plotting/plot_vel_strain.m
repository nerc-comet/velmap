%=================================================================
% plot_vel_strain.m
% Plot velocities and strains from velmap output.
% Tiledlayout requires 2019b or later.
%
% Modified from cahb_insarfit_gpsonly.m.
% Andrew Watson @ leeds, 12/07/2021
%=================================================================

%% setup

% main paths
resultpath = '/nfs/a285/homes/eearw/velmap_projects/out/khorrami_gnss_cleaned_smf-1.60/';
% gpspath = '/nfs/a285/homes/eearw/velmap_projects/gps/';
gpspath = '/scratch/eearw/decomp_frame_vels/gnss/khor/cleaned_stations/';

% files
gpsfile = 'khor_vert_10mm_gf7_buff01.csv';
velfile = 'velfit.dat';
strainfile = 'strain_savage_nring1.dat';
bordersfile = 'borderdata.mat';

% parameters
meshspacing = 0.05;
places = {'Iran Islamic Republic of','Iraq','Afghanistan','Turkey',...
    'Turkmenistan','Pakistan','Saudi Arabia','Armenia','Azerbaijan'}; % to plot

%% load inputs

% results
gps = readmatrix([gpspath gpsfile]);
vel = readmatrix([resultpath velfile]);
strain = readmatrix([resultpath strainfile]);

% plotting
borders = load('borderdata.mat');
load('vik.mat')
load('acton.mat'); acton = flipud(acton); % reverse colour progression
% plotting coords
lon = vel(:,1); lat = vel(:,2);

%% format velocities

[egrid,ngrid,ugrid] = format_vel(vel,meshspacing);

%% format strain

[strainmap,ssgrid,msgrid,i2grid] = format_strain(strain,meshspacing);

%% plot vels

clim = [-25 25];

figure()
tiledlayout(2,2,'TileSpacing','compact')

nexttile; hold on
plt_data(lon,lat,egrid,gps,strainmap.lonlims,strainmap.latlims,[-10 10],borders,places,'East vel (mm/yr)')
colormap(vik)

nexttile; hold on
plt_data(lon,lat,ngrid,gps,strainmap.lonlims,strainmap.latlims,clim,borders,places,'North vel (mm/yr)')
colormap(vik)

nexttile; hold on
plt_data(lon,lat,ugrid,gps,strainmap.lonlims,strainmap.latlims,clim,borders,places,'Up vel (mm/yr)')
colormap(vik)

%% plot strains

clim = [0 2e-7];

figure()
tiledlayout(2,2,'TileSpacing','compact')

nexttile; hold on
plt_data(lon,lat,i2grid,gps,strainmap.lonlims,strainmap.latlims,clim,borders,places,...
    'Second invariant of strain rate (/yr)')
colormap(acton)

nexttile; hold on
plt_data(lon,lat,msgrid,gps,strainmap.lonlims,strainmap.latlims,clim,borders,places,...
    'max shear strain rate (/yr)')
colormap(acton)

nexttile; hold on
plt_data(lon,lat,ssgrid,gps,strainmap.lonlims,strainmap.latlims,clim,borders,places,...
    'Savage and Simpson [1997] function of strain rate (/yr)')
colormap(acton)

%% plotting functions -----------------------------------------------------
function plt_data(lon,lat,data,gps,lonlim,latlim,clim,borders,places,titlestr)

% figure(); hold on

% plot input data
imagesc(lon,lat,data);

% add country borders
for ii = 1:length(places)
    b_ind = find(strcmp(borders.places,places(ii)));
    plot(borders.lon{b_ind},borders.lat{b_ind},'k')
end

% plot gps arrows
sc = 0.2; clat = 33.75; % scaling factor and centre lat
quiver(gps(:,1),gps(:,2),sc/cosd(clat)*gps(:,3),sc*gps(:,4),0,'color','k','linewidth',0.5);

axis xy
xlim(lonlim)
ylim(latlim)

colorbar
caxis(clim)

title(titlestr)

end
