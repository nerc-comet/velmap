%=================================================================
% difference_results.m
% Load two results directories and calculate the difference between the
% results.
% Tiledlayout requires 2019b or later.
%
% Andrew Watson @ leeds, 17/02/2022
%=================================================================

%% setup

% first results dir
resultpath1 = '/nfs/a285/homes/eearw/velmap_projects/out/khorrami_gnss_fine_smf-1.60/';
velfile1 = 'velfit.dat';
strainfile1 = 'strain_savage_nring1.dat';

% second results dir
resultpath2 = '/nfs/a285/homes/eearw/velmap_projects/out/zagros_gps_insar_025_smf-1.60/';
velfile2 = 'velfit.dat';
strainfile2 = 'strain_savage_nring1.dat';


% 
% % second results dir
% resultpath2 = '/nfs/a285/homes/eearw/velmap_projects/out/chris_gnss_fine_smf-1.60/';
% velfile2 = 'velfit.dat';
% strainfile2 = 'strain_savage_nring1.dat';



% % first results dir
% resultpath1 = '/nfs/a285/homes/eearw/velmap_projects/out/chris_gnss_fine_smf-1.60/';
% velfile1 = 'velfit.dat';
% strainfile1 = 'strain_savage_nring1.dat';

% second results dir
resultpath2 = '/nfs/a285/homes/eearw/velmap_projects/out/chris_cleaned_gnss_fine_smf-1.60/';
velfile2 = 'velfit.dat';
strainfile2 = 'strain_savage_nring1.dat';



% borders for plotting
bordersfile = 'borderdata.mat';
places = {'Iran Islamic Republic of','Iraq','Afghanistan','Turkey',...
    'Turkmenistan','Pakistan','Saudi Arabia','Armenia','Azerbaijan'}; % to plot

% mesh spacing
meshspacing = 0.05;

%% load files

% first results
vel1 = readmatrix([resultpath1 velfile1]);
strain1 = readmatrix([resultpath1 strainfile1]);

% second results
vel2 = readmatrix([resultpath2 velfile2]);
strain2 = readmatrix([resultpath2 strainfile2]);

% plotting
borders = load('borderdata.mat');
load('vik.mat')
load('acton.mat'); acton = flipud(acton); % reverse colour progression

% plotting coords
lon = vel1(:,1); lat = vel1(:,2);

%% format velocities and strains

% first
[egrid1,ngrid1,ugrid1] = format_vel(vel1,meshspacing);
[strainmap1,ssgrid1,msgrid1,i2grid1] = format_strain(strain1,meshspacing);

% second
[egrid2,ngrid2,ugrid2] = format_vel(vel2,meshspacing);
[strainmap2,ssgrid2,msgrid2,i2grid2] = format_strain(strain2,meshspacing);

%% calc diffs

egrid_diff = egrid1-egrid2;
ngrid_diff = ngrid1-ngrid2;
ugrid_diff = ugrid1-ugrid2;
ssgrid_diff = ssgrid1-ssgrid2;
msgrid_diff = msgrid1-msgrid2;
i2grid_diff = i2grid1-i2grid2;

%% plot vels

% east vel
clim = [-10 10];
plt_triplet(lon,lat,egrid1,egrid2,egrid_diff,strainmap1.lonlims,strainmap1.latlims,...
    clim,borders,places,'East vel (mm/yr)')
colormap(vik)

% north vel
clim = [-25 25];
plt_triplet(lon,lat,ngrid1,ngrid2,ngrid_diff,strainmap1.lonlims,strainmap1.latlims,...
    clim,borders,places,'North vel (mm/yr)')
colormap(vik)

%% plot strains

clim = [0 2e-7];

% 2nd inv
plt_triplet(lon,lat,i2grid1,i2grid2,i2grid_diff,strainmap1.lonlims,strainmap1.latlims,...
    clim,borders,places,'2nd invariant (/yr)')
colormap(acton)

% shear
plt_triplet(lon,lat,msgrid1,msgrid2,msgrid_diff,strainmap1.lonlims,strainmap1.latlims,...
    clim,borders,places,'Max shear (/yr)')
colormap(acton)

%% plotting functions -----------------------------------------------------
function plt_data(lon,lat,data,lonlim,latlim,clim,borders,places,titlestr)

% plot input data
imagesc(lon,lat,data);

% add country borders
for ii = 1:length(places)
    b_ind = find(strcmp(borders.places,places(ii)));
    plot(borders.lon{b_ind},borders.lat{b_ind},'k')
end

axis xy
xlim(lonlim)
ylim(latlim)

colorbar
caxis(clim)

title(titlestr)

end

function plt_triplet(lon,lat,data1,data2,data_diff,lonlim,latlim,clim,borders,places,titlestr)

figure()
tiledlayout(1,3,'TileSpacing','compact')

nexttile; hold on
plt_data(lon,lat,data1,lonlim,latlim,clim,borders,places,[titlestr ' - first'])

nexttile; hold on
plt_data(lon,lat,data2,lonlim,latlim,clim,borders,places,[titlestr ' - second'])

nexttile; hold on
plt_data(lon,lat,data_diff,lonlim,latlim,clim,borders,places,[titlestr ' - diff'])

end