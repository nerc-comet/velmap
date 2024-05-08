function [par,gpspar,insarpar,smpar,tssmpar,profpar] = readparfile(cfgfile,legacy_mode)
%=================================================================
% function [parmat] = readparfile(cfgfile)
%-----------------------------------------------------------------
% Function to get parameter matrix from a configure file
%                                                                  
% INPUT:                                                           
%   cfgfile: path to parameter text file (e.g. velmap.conf)
%   legacy_mode: toggles between old format (ninsarfile and insardir_n) and 
%                new format (insardir).
% OUTPUT:                                                          
%   par:  structure containing general parameters
%   gpspar: structure containing gps parameters
%   insarpar: structure containing insar parameters
%   smpar: structure containing spatial smoothing parameters
%   tssmpar: structure containing time series smoothing parameters
%   
% This merges the original readparfile.m, getpars_velmap.m, and getpar.m.
%
% Hua Wang @ Uni Leeds, 20/07/2008
% Andrew Watson @ Leeds, 15/06/2021
%                                                                  
% NOTE: use '#' for comments in config file, and ':' to seperate names and
% values (e.g. inv_e:   1)
%=================================================================

%% open config file

% test that config file exists
if ~isfile(cfgfile)
    disp('Config file not found, continuing with default values.')
end

% load the config file as a cell array
cfgcell = readcell(cfgfile,'FileType','text','CommentStyle','#');

%% format setup from cell array to structure

% processing mode
par.procmode = getparval(cfgcell,'procmode',1);

% output velocity parameters
par.invenu(1) = getparval(cfgcell,'inv_e',1);
par.invenu(2) = getparval(cfgcell,'inv_n',1);
par.invenu(3) = getparval(cfgcell,'inv_u',1);

% output directory
par.outdir = getparval(cfgcell,'outdir');

% mesh file
par.meshfile = getparval(cfgcell,'meshfile');

% strain
par.nring = getparval(cfgcell,'nring',1);

%% gps

% gps data
gpspar.ngpsfile = getparval(cfgcell,'ngpsfile',1);
for ii = 1:gpspar.ngpsfile
    gpspar.gpsfiles{ii} = getparval(cfgcell,['gpsfile_' num2str(ii)]);
    gpspar.gpsvel3d{ii} = getparval(cfgcell,['gpsvel3d_' num2str(ii)],0);
end

%% insar

% insar data
if legacy_mode == 1
    insarpar.ninsarfile = getparval(cfgcell,'ninsarfile');
    for ii = 1:insarpar.ninsarfile
        insarpar.dir{ii} = getparval(cfgcell,['insardir_' num2str(ii)]);
    end
else
    insarpar.ninsarfile = sum(strcmp(cfgcell(:,1),'insardir'));
    for ii = 1:insarpar.ninsarfile
        insarpar.dir{ii} = getparval(cfgcell,'insardir',[],ii);
    end
end

% insar and errors file extensions
insarpar.insar_ext = getparval(cfgcell,'insar_ext','.vel.geo.tif');
insarpar.errors_ext = getparval(cfgcell,'insar_ext','.vstd.geo.tif');

% insar resolution
insarpar.xpsize = getparval(cfgcell,'insar_xpsize',0);
insarpar.ypsize = getparval(cfgcell,'insar_ypsize',0);

if insarpar.xpsize == 0 || insarpar.ypsize == 0
  insarpar.lksx = getparval(cfgcell,'insar_lksx',10);
  insarpar.lksy = getparval(cfgcell,'insar_lksy',10);
end

% orbital fitting parameter
insarpar.orbdegree = getparval(cfgcell,'orbdegree',1);

% atm fitting parameter
insarpar.atmdegree = getparval(cfgcell,'atmdegree',1);

% inversion components
insarpar.invenu = [getparval(cfgcell,'inv_e',1)...
                    getparval(cfgcell,'inv_n',1)...
                     getparval(cfgcell,'inv_u',1)];

%% spatial smoothing
% smfactor: smoothing factor(0: calculate & plot L-curve; others: smoothing factor)
% smf_min/max/int: region of smoothing factors for L-curve calculation, 
% the exact region will be calculated by 10^(smf)
% lcurve_lksx/lksy: looks number for L-curve calculation

smpar.smf = getparval(cfgcell,'smfactor');

if smpar.smf == 0 || smpar.smf == 999
  smpar.smf_min = getparval(cfgcell,'smf_min');
  smpar.smf_max = getparval(cfgcell,'smf_max');
  smpar.smf_int = getparval(cfgcell,'smf_int');
  smpar.lcurv_lksx = getparval(cfgcell,'lcurv_lksx');
  smpar.lcurv_lksy = getparval(cfgcell,'lcurv_lksy');
  
  if (smpar.smf_min > smpar.smf_max)
    error('smoothing factors are wrong, check smf_min/smf_max');
  end
  
else
  smpar.smf=10^smpar.smf;
  
end

%% temporal domain smoothing parameters
% smfactor: smoothing factor(0: calculate & plot L-curve; others: smoothing factor)
% smf_min/max/int: region of smoothing factors for L-curve calculation, 
% the exact region will be calculated by 10^(smf)
% lcurve_lksx/lksy: looks number for L-curve calculation
% mingps: mininum number of gps sites for each epoch

if par.procmode == 2
  tssmpar.t0 = getparval(cfgcell,'tst0');   %yyyymmdd
  tssmpar.t0 = datenum(num2str(tssmpar.t0),'yyyymmdd'); %serial days
  tssmpar.dt = getparval(cfgcell,'tsdt');   %in days
  tssmpar.smf = getparval(cfgcell,'tssmfactor'); 
  tssmpar.smorder = getparval(cfgcell,'tssmorder',2); 
  
  if tssmpar.smf==0 || tssmpar.smf==999
    tssmpar.smf_min = getparval(cfgcell,'tssmf_min');
    tssmpar.smf_max = getparval(cfgcell,'tssmf_max');
    tssmpar.smf_int = getparval(cfgcell,'tssmf_int');
    
    if (tssmpar.smf_min > tssmpar.smf_max)
      error('smoothing factors are wrong, check tssmf_min/tssmf_max');      
    end
    
  else
    tssmpar.smf=10^tssmpar.smf;
    
  end
  
  tssmpar.mingps = getparval(cfgcell,'tsmingps',10);
  
else    
    tssmpar = nan;
  
end

%% profiles

profpar.profflag = getparval(cfgcell,'make_prof');

if profpar.profflag==1
  profpar.swath = getparval(cfgcell,'profswath');
  profpar.step = getparval(cfgcell,'profstep');
  profpar.gpsswath = getparval(cfgcell,'profswath_gps',50);
  profpar.insarswath = getparval(cfgcell,'profswath_insar',50);
  
else
    profpar = nan;
    
end

%make velocity field on a regular grid
% grdvel.dx = getparval(cfgcell,'grdveldx');
% grdvel.dy = getparval(cfgcell,'grdveldy');

%% getparval ==============================================================
% cfgcell: n-by-2 cell array containing parameters
% parstring: name of desired parameter
% defval: default val if parameter is not found
% indn: index to use when multiple matches are found
function [val] = getparval(cfgcell,parstring,defval,indn)
    
    if nargin ~= 4
        indn = 1;
    end
    
    try
%         val = cfgcell{strcmp(cfgcell(:,1), parstring),2};
        indx = find(strcmp(cfgcell(:,1), parstring));
        val = cfgcell{indx(indn),2};
    catch
        val = defval;
    end
end

end