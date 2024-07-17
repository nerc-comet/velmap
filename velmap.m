%=================================================================
% velmap
% Main script to make a crustal velocity map using geodetic data.
%
% Calls functions from the main, pilib, and plotting directories.
% Written for Matlab 2019a.
%
% Can use the original parameter file setup or a new version. Toggle for
% this is found in readparfile.
%
% Original - Hua Wang @ Leeds, 17/09/2009
% Updated - Andrew Watson @ leeds, 15/06/2021
% Modified - Dehua Wang @ leeds, 03/04/2023
%=================================================================

%% 0. setup

% input config file
cfgfile = 'examples/ATF/ATF_gps_insar.conf';

% add path to functions
addpath('main')

%% 1. read config file

fprintf('===> reading config file ...\n');
% 1 = legacy mode, 0 = new parameter file setup
[par,gpspar,insarpar,smpar,tssmpar,profpar] = readparfile(cfgfile,0);

%% 2. make triangular mesh
% fault location, range of study area needed

fprintf('===> reading mesh ...\n');

if exist(par.meshfile,'file')
    if strcmp(par.meshfile(end-3:end),'.msh')
        trim=gid2mat(par.meshfile);
        
    elseif strcmp(par.meshfile(end-3:end),'.mat')
        load(par.meshfile)
        
    end
else
    error('mesh file is not available');
    
end

%% 3. prepare gps observations

fprintf('===> loading gps data ...\n');
if gpspar.ngpsfile>0    
    [gps]=loadgps(gpspar);        %load gps data
    
    for i=1:gpspar.ngpsfile
        [gps(i).site]=tidygps(trim,gps(i).site); %remove gps outside of the mesh
        gps(i).nsite=length(gps(i).site);
    end
    
else
    fprintf('===> no gps data to load...exiting...\n');
    return;
    
end

%% 4. prepare insar observations

if insarpar.ninsarfile>0
    fprintf('===> loading insar data ...\n');
%  [insar]=loadinsar(insarpar);
    [insar]=loadlics(insarpar);  %load insar data in lics format
  
else
    fprintf('===> no insar data to load...\n');
    insar=[];
    
end

plottrim(trim,gps,insar);
drawnow;


%% 5. find best smoothing factor

if smpar.smf==0
  fprintf('===>  processing for all smoothing factors...\n');
  smpar.smf = (smpar.smf_min:smpar.smf_int:smpar.smf_max);
  smpar.smf = 10.^smpar.smf;
  
elseif smpar.smf==999 % currently non-functioning
  fprintf('===>  find the best smoothing factor...\n');
  smpar.smf = lcurve_vmp(trim,smpar,gps,insar,1,par.outdir);
  
end


%% solution for each smoothing factor

%update invenu
% invenu=getinvenu(gps,insar);
invenu = insarpar.invenu;
nsmf=length(smpar.smf);

for i=1:nsmf
  
  fprintf('===> processing smoothing factor %d/%d\n',i,nsmf);
  %output directory
  smfdir=char(strcat(par.outdir,'smf',num2str(log10(smpar.smf(i)),'%+4.2f'),'/'));
  if ~exist(smfdir,'dir')
    mkdir(smfdir)
  end
  
  %% 6. solve system of equations
  
  [fitmodel,vcmmodel,wrss,rough] = ...
      solve_vmp(trim,smpar.smf(i),gps,insar,smfdir,invenu,0);
  lcv = [smpar.smf(i) rough wrss log10(smpar.smf(i))]; %JRW add
  
  lcvfile=strcat(smfdir,'lcv.txt');
  dlmwrite(lcvfile,lcv,'precision','%4.4f','delimiter','\t')
  
  % JRW add
  writematrix(lcv,[smfdir 'lcv.dat'],'delimiter','\t')
  save([smfdir 'fitmodel'],'fitmodel','-v7.3');
  save([smfdir 'vcmmodel'],'vcmmodel','-v7.3');
  

  %% 7. forward calculation

  % output fitted velocity field
  fprintf('===> output fitted velocity field... \n');
  fitvtx = fitmodel2vel(trim,fitmodel,vcmmodel,invenu,smfdir);
  
  save([smfdir 'fitvtx'],'fitvtx','-v7.3');

  % forward calculation for gps
  fprintf('===> forward calculating... \n');
  gpsfit = gpsfwd(trim,fitmodel,vcmmodel,invenu,gps,smfdir);
  
  save([smfdir 'gpsfit'],'gpsfit','-v7.3');
  
  % forward calculation for insar
  if insarpar.ninsarfile>0
    insarfit = insarfwd(insar,trim,fitmodel,invenu,smfdir,gps);
    save([smfdir 'insarfit'],'insarfit','-v7.3');    
  end

  save precrash.mat

  %% 8. strain rate

  % calculate strain rate
  fprintf('===> calculating strain rate... \n');
  nvtx=length(trim.x);
  fitvel=(reshape([fitvtx.vel],[],nvtx))';
  strain = vel2strain_tri(trim,fitvel,smfdir);
  [strain,eulervec] = vel2strain_savage(trim,fitvel,vcmmodel,par.nring,smfdir,2);
% 
%   %-------------------
%   %9. profile
%   %-------------------
%   %make profile for the velocity field
%   if profflag==1
%     fprintf('===> making profiles... \n');
%     %interactively extract profile
%     profdef=char('profdef.dat');
%     if ~exist(profdef,'file')
%       plotvel(fitvtx,gps,gpsfit,faults,smfdir);
%       extractprof(prof.swath,prof.step);
%     end
% 
%     %read profile
%     [prof]=profdefine(profdef);
% 
%     %extract fault position on the profile
%     if ~exist('faultonprof.dat','file')
%       extractfaultonprof(prof,faults);
%     end
% 
%     %calculate profile
%     nprof=length(prof);
%     profdir=strcat(smfdir,'prof/');
%     if ~exist(profdir,'dir')
%       mkdir(profdir);
%     end
%     for iprof=1:nprof
%       fprintf('making profile %d/%d\n',iprof,nprof);
%       if prof(iprof).swath==0
%         make_profline_vel(trim,fitmodel,vcmmodel,invenu,prof(iprof),profdir);
%       else
%         make_profswath_vel(trim.x,trim.y,fitvel,vcmmodel(1:2*nvtx,1:2*nvtx),prof(iprof),profdir);
%       end
%       %make profile for the observed GPS data
%       %low efficiency to make profile for each site once a time
%       gpsprofdef=prof(iprof);
%       gpsprof=[];
%       gpsprofdef.swath=gpsswath;
%       for igf=1:gpspar.ngpsfile
%         ns=length(gps(i).site);
%         for is=1:ns
%           isite=gps(igf).site(is);
%           igpsprof=make_profswath_vel(isite.lon,isite.lat,isite.vel,isite.vcm,gpsprofdef);
%           gpsprof=[gpsprof;igpsprof];
%         end
%       end
%       if size(gpsprof,1)>0
%         outfile=strcat(profdir,gpsprofdef.id,'.prof_gps');
%         save(char(outfile),'gpsprof','-ASCII');
%       end
% 
%       %make profile for the observed InSAR data
%       insarprofdef=prof(iprof);
%       insarprofdef.swath=insarswath;
%       for isf=1:insarpar.ninsarfile
%         %stackmap=insar(isf).stackmap; 
%         %using original resolution stackmap
%         %ifghdr=rsc2hdr(char(strcat(insarpar.dir(isf),'ratemap/ifg.rsc')));
%         %ifghdr=rsc2hdr(char(strcat(insarpar.dir(isf),'ifghdr.mat')));
%         ifghdr=char(strcat(insarpar.dir(isf),'ifghdr.mat')); %old version of pi-rate 2.0 or earlier - JRW add
%         load(ifghdr); %JRW add
%         %stackmap=readmat(char(strcat(insarpar.dir(isf),'ratemap/stackmap.dat')),ifghdr.length,ifghdr.width,1);
%         stackmap=readmat(char(strcat(insarpar.dir(isf),'stackmap.dat')),ifghdr.length,ifghdr.width,1); %JRW add
%         [sarprof_pt,sarprof] = make_prof(stackmap,insarprofdef,ifghdr);
%         if size(sarprof,1)>0
%           sarprof_pt=double(sarprof_pt);
%           sarprof=double(sarprof);
%           outfile=strcat(profdir,num2str(insarprofdef.id),'.prof_insar',num2str(isf,'%02d'));
%           save(char(outfile),'sarprof','-ASCII');
%           outfile=strcat(profdir,num2str(insarprofdef.id),'.prof_insar',num2str(isf,'%02d'),'_pt');
%           save(char(outfile),'sarprof_pt','-ASCII');
%         end
%       end
%     end
%   end
% 
% %make velocity field on a regular grid
% disp('===> NOT making velocity field on a regular grid...')
% %make_grid_vel(trim,fitmodel,vcmmodel,invenu,grdvel.dx,grdvel.dy,smfdir);
end
% %remove pars.mat
% !rm -f pars.mat
% fprintf('====finished successfully, congratulations!====\n');
