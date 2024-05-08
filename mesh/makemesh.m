%-------------------------------------------------
% makemesh
%
% Updated version of the original makemesh.m.
% Modified to be a seperate script ran before velmap, as opposed to a
% function.
% Includes low res zones.
% 
% Hua Wang, 29/09/2015
% Disable saved, added randomisation - Andrew Watson, 19/07/2021
%-------------------------------------------------

%% toggles

lowreszones = 0; % (0,1) read polygons for lower resolution areas
meshoutline = 0; % (0,1) read polygon for mesh outline

preview = 1; % (0,1) simple mesh plot
save_mesh = 1; % (0,1) save mesh to .mat file

%% setup

% in files
polyfile = '/nfs/a285/homes/eearw/velmap_projects/mesh/lowrespoly.txt';
outlinefile = '/nfs/a285/homes/eearw/velmap_projects/mesh/outline.txt';

% out files
outfile = '/nfs/a285/homes/eearw/velmap_projects/mesh/net_iran_simple_025.mat';

% x limits and intervals for full and low res mesh
xmin = 42;
xmax = 67;
dx = 0.25;
dxlow = 0.5;

% y limits and interval for full and low res mesh
ymin = 24;
ymax = 42.5;
dy = 0.25;
dylow = 0.5;

% maximum fractions of a degree added or subtracted from locations to randomize
dr = dx/5;

%% generate points

% generate uniform grids
[x,y]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);

x = x(:); y = y(:);

%% add low res zone

if lowreszones == 1
    
    % load polygons
    poly = readpoly(polyfile);
    
    % define low res grid
    [xlow,ylow] = meshgrid(xmin:dxlow:xmax,ymin:dylow:ymax);
    xlow = xlow(:); ylow = ylow(:);
    
    % for each polygon...
    for ii = 1:length(poly)
        
        % check high res points inside polygon
        [in_poly,~] = inpolygon(x,y,poly{ii}(:,1),poly{ii}(:,2));
        
        % remove those
        x(in_poly) = []; y(in_poly) = [];
        
        % check low res points
        [in_poly,~] = inpolygon(xlow,ylow,poly{ii}(:,1),poly{ii}(:,2));
        
        % add onto high res points
        x = [x; xlow(in_poly)]; y = [y; ylow(in_poly)];
    end
    
end

%% add random shift

% add random shift scaled from 0 to dr
ry = dr .* rand(size(y));
rx = dr .* rand(size(x));

% ignore edge points
bind = boundary(x,y);
rx(bind) = 0; ry(bind) = 0;

x = x + rx; y = y + ry;

%% mask points outside of outline

if meshoutline == 1
    
    % load polygon
    outline = readpoly(outlinefile);
    
    % check points inside outline
    [in_poly,~] = inpolygon(x,y,outline(:,1),outline(:,2));
    
    % remove those outside from mesh
    x(~in_poly) = []; y(~in_poly) = [];
    
    % get boundary inds and conver to edge list
    bind = boundary(x,y);
    bind_shift = circshift(bind,-1);
    outline = [bind(1:end-1) bind_shift(1:end-1)];
%     outline = reshape(boundary(x,y),2,[])';
    
end

%% generate triangles

% addpath(genpath('/nfs/a285/homes/eearw/mesh2d'))
% [vert,etri,tria,tnum] = refine2([x y]);%,outline) ;
% 
% trim.tri = tria;
% trim.x = vert(:,1);
% trim.y = vert(:,2);

% calculate delaunary triangles
% trim.tri = delaunayTriangulation(x,y,outline);
trim.tri = delaunay(x,y);
trim.y=reshape(y,[],1);
trim.x=reshape(x,[],1);

%% mask triangles outside of outline

% if meshoutline == 1
%     
%     % load polygon
%     outline = readpoly(outlinefile);
%     
%     % check points inside outline
%     [in_outline,~] = inpolygon(trim.x,trim.y,outline(:,1),outline(:,2));
%     
%     % get indices
%     outline_ind = find(~in_outline);
%     
%     % remove points outside outline
%     trim.x(~in_outline) = []; trim.y(~in_outline) = [];
%     
%     % remove triangle with any point outside outline
%     trim.tri(any(ismember(trim.tri,outline_ind),2), :) = [];
%     
% end

%% preview mesh

if preview == 1
    figure()
    triplot(trim.tri,trim.x,trim.y,'color',[0.625 0.625 0.625],'linewidth',0.25);
end

%% output

if save_mesh == 1
    save(outfile,'trim')
end
