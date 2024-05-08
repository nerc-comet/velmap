function [] = plottrim(trim,gps,insar)
%===========================================
%function [] = plottrim(trim,gps,insar)
%
% Plot triangular mesh and insar data together
%
% Input:
%  trim: triangular mesh
%  gps:  gps data
%  insar: insar data
%
% Output:
%  None
%
% Hua Wang @ Uni Leeds, 01/10/2009
%
% 23/09/2015 HW: revise for new gps structure
%                renamed as plotrim as plotmesh is a embeded function
%===========================================

%plot triangular mesh
figure

if nargin>2
  ninsar=length(insar);
  for is=1:ninsar
    [rows,cols]=size(insar(is).stackmap);
    n=rows*cols;
    vrate=reshape(insar(is).stackmap',n,1);
    [xx,yy]=meshgrid(1:cols,1:rows);
    xxv=reshape(xx',n,1);
    yyv=reshape(yy',n,1);
    xxv = insar(is).ifghdr.xfirst+(xxv-1)*insar(is).ifghdr.xstep;
    yyv = insar(is).ifghdr.yfirst+(yyv-1)*insar(is).ifghdr.ystep;
    clear('xx','yy');
    xxv(isnan(vrate))=[];
    yyv(isnan(vrate))=[];
    vrate(isnan(vrate))=[];
    r=is/ninsar;
    g=1-is/ninsar;
    b=1-is/ninsar;
    plot(xxv,yyv,'.','MarkerSize',10,'MarkerEdgeColor',[r g b]);
    hold on
  end
end

triplot(trim.tri,trim.x,trim.y);

hold on

if nargin>1
  if isfield(gps,'site')
    site=gps.site;
  else
    site=gps;
  end

  scale=50;
  ngps=length(site);
  vel=[];
  for i=1:ngps
    if isfield(site(i),'vel')
      vel=[vel;site(i).vel];
    else
      vel=[vel;site(i).tsvel(1,:)];
    end
  end
  quiver([site.lon]',[site.lat]',vel(:,1)/scale,vel(:,2)/scale,'k','AutoScale','off');
end

axis equal
axis([min(trim.x)-1,max(trim.x)+1,min(trim.y)-1,max(trim.y)+1]);

