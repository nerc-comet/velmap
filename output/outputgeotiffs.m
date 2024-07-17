%=================================================================
% outputgeotiffs.m
%
%script to output geotiffs for each frame
%
%
% Originally from Tim Wright
% Modified by Dehua Wang, University of Leeds, 04/03/2023
%=================================================================

load examples/ATF/output/smf-1.10/precrash.mat
%% read in files produced by insarfit and the original interferograms
load examples/ATF/output/smf-1.10/insar_1lk.mat
load examples/ATF/output/smf-1.10/insarfit2.mat
mkdir geotiffs-1.1-m
ninsar=length(insar);

for i=1:ninsar
   
    % read original velocity file to get header
    namestruct=dir(string(strcat(insarpar.dir(i),'vstd.geo.tif')));
    stackmapname=sprintf("%s/%s", string(insarpar.dir(i)),namestruct.name);
    [vstd,R]=geotiffread(stackmapname);
    
    %write velocities referenced to Eurasia
    filename = sprintf("%s%d%s","geotiffs-1.1-m/vel_eurasiaref_frame_",i,".tif");
    outgrid = insarfit2(i).ratemap+insarfit2(i).resmap; %this is the original interferogram with atmosphere and orbital correction terms applied to put it into geocoded coordinates.
%   outgrid = insar(i).stackmap-insarfit2(i).orbmap;
    %outgrid(isnan(outgrid))= -9999;
    outgrid = -outgrid;      % reverse the direction, Dehua Wang
    geotiffwrite(filename,outgrid,R);

    %write vstd
    filename = sprintf("%s%d%s","geotiffs-1.1-m/vstd_",i,"vstd.geo.tif");
    %vstd(isnan(vstd))= -9999;
    geotiffwrite(filename,vstd,R);  %just rewriting what we read in to make life easier later on

    %write orbmap
    filename = sprintf("%s%d%s","geotiffs-1.1-m/orbmap_",i,".tif");
    outgrid = insarfit2(i).orbmap;
    %outgrid(isnan(outgrid))= -9999;
    geotiffwrite(filename,outgrid,R);
    
    %write atmmap
    filename = sprintf("%s%d%s","geotiffs-1.1-m/atmmap_",i,".tif");
    outgrid = insarfit2(i).atmmap;
    %outgrid(isnan(outgrid))= -9999;
    geotiffwrite(filename,outgrid,R);
    
    %write resmap
    filename = sprintf("%s%d%s","geotiffs-1.1-m/resmap_",i,".tif");
    outgrid = insarfit2(i).resmap;
    outgrid = -outgrid;
    %outgrid(isnan(outgrid))= -9999;
    geotiffwrite(filename,outgrid,R);
    
    
    %write modeled los rate
    filename = sprintf("%s%d%s","geotiffs-1.1-m/modellosmap_",i,".tif");
    outgrid = insarfit2(i).ratemap;
    outgrid = -outgrid;
    %outgrid(isnan(outgrid))= -9999;
    geotiffwrite(filename,outgrid,R);

    
end
