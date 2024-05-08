%=================================================================
% plot_mesh.m
%
% Plots the mesh, gps velocities, and outlines of insar frames.
%
% Tiledlayout requires 2019b or later.
%
% Modified from plot_mesh.m.
% Andrew Watson @ leeds, 19/07/2021
%=================================================================

%% setup

% plotting toggles (1 = on, 0 = off)
plt_insar = 0;
plt_gnss = 0;

% script paths
addpath ../main

% main paths
resultpath = '/nfs/a285/homes/eearw/velmap_projects/out/zagros_gps_insar_smf-1.40/';
gpspath = '/nfs/a285/homes/eearw/velmap_projects/gps/';
meshpath = '/nfs/a285/homes/eearw/velmap_projects/mesh/';

% files
insarfitfile = 'insarfit.mat';
gpsfile = 'iran_tol1.5_minocc2.5_dist5_2D.dat';
meshfile = 'net_iran3.mat';
bordersfile = 'borderdata.mat';

% parameters
lonlim = [41 68.5];
latlim = [23 43.5];
sc = 0.1; clat = 33.75; % scaling factor and centre lat for gnss
places = {'Iran Islamic Republic of','Iraq','Afghanistan','Turkey',...
    'Turkmenistan','Pakistan','Saudi Arabia','Armenia','Azerbaijan'}; % to plot

%% load inputs

% data
if strcmp(meshfile(end-3:end),'.msh')
    trim = gid2mat([meshpath meshfile]);
else
    load([meshpath meshfile])
end

if plt_insar == 1; load([resultpath insarfitfile]); end
if plt_gnss == 1; gps = readmatrix([gpspath gpsfile]); end

% plotting
borders = load('borderdata.mat');

%% get insar boundaries

if plt_insar == 1
    
    insarboundary = cell(1,length(insarfit));

    for ii = 1:length(insarfit)

        % lon and lat coord vectors
        x = insarfit(ii).ifghdr.xfirst ...
            + [0:(insarfit(ii).ifghdr.width-1)].*insarfit(ii).ifghdr.xstep;
        y = insarfit(ii).ifghdr.yfirst ...
            + [0:(insarfit(ii).ifghdr.length-1)].*insarfit(ii).ifghdr.ystep;

        % boundary of each frame
        [xx,yy] = meshgrid(x,y);
        xx = xx(:); yy = yy(:);
        xx(isnan(insarfit(ii).ratemap(:))) = []; yy(isnan(insarfit(ii).ratemap(:))) = [];

        boundaryind = boundary(xx(:),yy(:));
        insarboundary{ii} = [xx(boundaryind) yy(boundaryind)];

    end

end

%% plot

figure(); hold on

% add country borders
for ii = 1:length(places)
    b_ind = find(strcmp(borders.places,places(ii)));
    plot(borders.lon{b_ind},borders.lat{b_ind},'r','LineWidth',2)
end

% mesh
triplot(trim.tri,trim.x,trim.y,'color',[0.625 0.625 0.625],'linewidth',0.25);

% gnss
if plt_gnss == 1
    quiver(gps(:,1),gps(:,2),sc/cosd(clat)*gps(:,3),sc*gps(:,4),0,'color','k','linewidth',0.5);
end

% insar boundaries
if plt_insar == 1
    for jj = 1:length(insarboundary)
        plot(insarboundary{jj}(:,1),insarboundary{jj}(:,2),'b')
    end
end

xlim(lonlim)
ylim(latlim)
