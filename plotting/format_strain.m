function [strainmap,ssgrid,msgrid,i2grid] = format_strain(strain,meshspacing)
%=================================================================
% format_strain.m
% Provide matrix of strains. Calculates components and returns regular
% grids for plotting.
%
% Andrew Watson @ leeds, 17/02/2022
%=================================================================

% extract components
strainmap.lon = strain(:,1);
strainmap.lat = strain(:,2);
strainmap.exx = strain(:,3);
strainmap.eyy = strain(:,4);
strainmap.exy = strain(:,5);
strainmap.dilat = strain(:,6);
strainmap.maxshear = strain(:,7);
strainmap.I2 = strain(:,11);

% mesh
[longrid,latgrid] = meshgrid(min(strainmap.lon):meshspacing:max(strainmap.lon),...
    min(strainmap.lat):meshspacing:max(strainmap.lat));

% calculate shear strain and dilatancy
strainmap.shear = sqrt((strainmap.exx - strainmap.eyy).^2 + 4*strainmap.exy.^2);
strainmap.dil = strainmap.exx + strainmap.eyy;

strainmap.eps1 = (strainmap.dil + strainmap.shear)/2;
strainmap.eps2 = (strainmap.dil - strainmap.shear)/2;

strainmap.savagesimpson ...
    = max([abs(strainmap.eps1),abs(strainmap.eps2),abs(strainmap.eps1 + strainmap.eps2)],[],2);

% interpolate onto uniform grid
ssinterp = scatteredInterpolant(strainmap.lon,strainmap.lat,strainmap.savagesimpson,'linear','none');
ssgrid = ssinterp(longrid,latgrid);

msinterp = scatteredInterpolant(strainmap.lon,strainmap.lat,strainmap.maxshear,'linear','none');
msgrid = msinterp(longrid,latgrid);

i2interp = scatteredInterpolant(strainmap.lon,strainmap.lat,strainmap.I2,'linear','none');
i2grid = i2interp(longrid,latgrid);

% calculate boundary of strain points
k = boundary(strainmap.lon,strainmap.lat,1);
gridinboundary = inpolygon(longrid,latgrid,strainmap.lon(k),strainmap.lat(k));

% mask points outside of boundary
ssgrid(gridinboundary==0) = NaN;
msgrid(gridinboundary==0) = NaN;
i2grid(gridinboundary==0) = NaN;

% strainmap.lonlims = [floor(min(strainmap.lon)) ceil(max(strainmap.lon))];
% strainmap.latlims = [floor(min(strainmap.lat)) ceil(max(strainmap.lat))];
strainmap.lonlims = [min(strainmap.lon) max(strainmap.lon)];
strainmap.latlims = [min(strainmap.lat) max(strainmap.lat)];

end

