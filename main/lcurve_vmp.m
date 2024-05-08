function [smfbest,wrss,rough,smf]=lcurve_vmp(trim,smpar,gps,insar,iplot,outdir)
%===============================================================
%function [smfbest,wrss,rough,smf]=lcurve_vmp(trim,smpar,gps,insar,iplot,outdir)
%                                                                    
% Calculate L-curve and best-fitted smoothing factor
%
% INPUT:
%  trim:      triangular mesh structure (tri,node)
%  smpar:     smoothing parameters
%  gps:       gps data
%  insar:     insar data structure
%  iplot:     flag to plot (1: plot, 0: no, default 1)
%  outdir:    output directory
%
% OUTPUT:
%  smfbest:   best-fitted smoothing factor
%  wrss:      weighted residual square sum
%  rough:     roughness
%  smf:       smoothing factors
%
% Hua Wang @ Uni Leeds, 28/09/2009
%
% 17/03/2015 HW: change filename to avoid conflict with pirate
%===============================================================

if nargin<4
  insar=[];
end
if nargin<5
  iplot=1;
end
if nargin<6
  outdir=[];
end

smf=(smpar.smf_min:smpar.smf_int:smpar.smf_max);
smf=(-2.2:.2:-0.8);
smf=10.^smf;
nsmf=length(smf);

%subsample the insar data on a sparse regular grid
ninsar=length(insar);
for i=1:ninsar
  ifghdr=insar(i).ifghdr;

  [insar(i).los,ix,iy]=downsmp(insar(i).los,smpar.lcurv_lksx,smpar.lcurv_lksy);
  vstackmap=reshape(insar(i).stackmap',numel(insar(i).stackmap),1);
  flag=(1:numel(insar(i).stackmap));
  flag(isnan(vstackmap))=[];
  idx=repmat(ix,length(iy),1)+repmat(((iy-1)*ifghdr.width)',1,length(ix));
  idx=reshape(idx',numel(idx),1);
  sel=[];
  for j=1:length(idx)
    isel = find(flag==idx(j));
    if ~isempty(isel)
      sel=[sel;isel];
    end
  end
  insar(i).vcm=insar(i).vcm(sel,sel');
  insar(i).azi=insar(i).azi(iy,ix);
  insar(i).stackmap=insar(i).stackmap(iy,ix);
  if insar(i).proc.atmdegree~=0
    insar(i).dem=insar(i).dem(iy,ix);
  end

  insar(i).ifghdr=ifghdrlooks(insar(i).ifghdr,smpar.lcurv_lksx,smpar.lcurv_lksy);
end

%calculate wrss and roughness
for i=1:nsmf
  fprintf('processing for the %d/%-d smoothing factor...\n',i,nsmf);
  [vel,stdvel,wrss(i),rough(i)]=solve_vmp(trim,smf(i),gps,insar,outdir,0);
end

%find best-fitted smoothing factor
%[smfbest]=lcorner(rough,wrss,smf,5,iplot);

%wrss .vs. roughness plot
if iplot==1
  figure
  plot(rough,wrss,'-bo','LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',6);
  xoff=(max(rough)-min(rough))*0.01;
  yoff=(max(wrss)-min(wrss))*0.02;
  for i=1:nsmf
    text(rough(i)+xoff,wrss(i)+yoff,num2str(log10(smf(i)),'%.1f'));
  end
  ylabel('weighted rms (mm)','fontsize',12)
  xlabel('solution roughness (mm/deg^2)','fontsize',12)
  %ylabel('weighted rss (mm^2)','fontsize',12)
  %xlabel('solution roughness (mm^2/deg^4)','fontsize',12)
  smfbest = input('PLEASE INPUT SMOOTHING FACTOR:');
  smfbest = 10^smfbest;
end

if ~isempty(outdir)
  %lcv=[smf' rough' wrss' log10(smf')];
  %save(char(strcat(outdir,'lcurve.dat')),'lcv','-ASCII');
  lcvfile=char(strcat(outdir,'lcurve.dat'));
  fid=fopen(lcvfile,'w');
  for i=1:nsmf
    fprintf(fid,'%10.4f %10.4f %10.4f %10.2f\n',smf(i),rough(i),wrss(i),log10(smf(i)));
  end
  fclose(fid);
end
