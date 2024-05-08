function [gps] = loadgps(gpspar)
%=============================================
%function [gps] = loadgps(gpspar)
%
% Load gps data from an ascii file
%
% Input:
%   gpspar: gps parameters
%
% Output:
%   gps=struct('ll',{},'vel',{},'vcm',{},'staid',{});
%
% Hua Wang @ Uni Leeds, 22/09/2009
% Andrew Watson @ Leeds, 16/06/2021
%
% 24/02/2022 AW: Add alternative file read incase of csv
% 16/06/2021 AW: updated to new gpspar structure
% 16/03/2015 HW: load multiple gps files
% 19/08/2011 HW: change the structure for time series
% 27/03/2010 HW: support 3D velocities
%=============================================
%read text file

for igf = 1:gpspar.ngpsfile
  %filename=char(gpspar.filename(igf));
  filename = gpspar.gpsfiles{igf};
  
  if gpspar.gpsvel3d{igf}==0
    gps(igf).invenu=[1 1 0];
    try % original method of opening gnss file
        [lon,lat,ve,vn,stde,stdn,coven,staid] ...
            = textread(filename,'%f %f %f %f %f %f %f %s','commentstyle','shell');
    catch % if that fails, assume csv input and try readmatrix
        gps_data = readmatrix(filename);
        
        % unpack into variables
        lon = gps_data(:,1); lat = gps_data(:,2);
        ve = gps_data(:,3); vn = gps_data(:,4);
        stde = gps_data(:,5); stdn = gps_data(:,6);
        
        % if no covariance is given, assume its zero
        if size(gps_data,2) >= 7
            coven = gps_data(:,7);
        else
            coven = zeros(size(gps_data,1),1); 
        end
        
        % read matrix will fail on station names, so just assume nans
        staid = nan(size(gps_data,1),1); 

        clear gps_data
    end
            
    
  else
    gps(igf).invenu=[1 1 1];
    try
        [lon,lat,ve,vn,vu,stde,stdn,stdu,coven,coveu,covnu,staid] ...
            = textread(filename,'%f %f %f %f %f %f %f %f %f %f %f %s','commentstyle','shell');
    catch
        % as above, but for 3D gnss
        gps_data = readmatrix(filename);
        lon = gps_data(:,1); lat = gps_data(:,2);
        ve = gps_data(:,3); vn = gps_data(:,4); vu = gps_data(:,5);
        stde = gps_data(:,6); stdn = gps_data(:,7); stdu = gps_data(:,8);
        coven = gps_data(:,9); coveu = gps_data(:,10); covnu = gps_data(:,11);
        staid = nan(size(gps_data,1),1); 
        clear gps_data
    end
    
  end

  %EN component
  gps(igf).ndim=sum(gps(igf).invenu);
  gps(igf).nsite=length(ve);
  
  for i=1:gps(igf).nsite
    gps(igf).site(i).lon=lon(i);
    gps(igf).site(i).lat=lat(i);
    gps(igf).site(i).staid=staid(i);
    coven(i)=(stde(i).*stdn(i)).*coven(i);
    
    if gpspar.gpsvel3d{igf}==0
      gps(igf).site(i).vel=[ve(i) vn(i)];
      gps(igf).site(i).vcm=[stde(i).^2, coven(i); coven(i), stdn(i).^2];
      
    else
      coveu(i)=(stde(i).*stdu(i)).*coveu(i);
      covnu(i)=(stdn(i).*stdu(i)).*covnu(i);
      gps(igf).site(i).vel=[ve(i), vn(i), vu(i)];
      gps(igf).site(i).vcm=[stde(i).^2, coven(i), coveu(i); coven(i), stdn(i).^2, covnu(i); coveu(i), covnu(i), stdu(i).^2];
      
    end
  end
end
