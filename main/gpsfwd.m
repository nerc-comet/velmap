function[gpsfit]=gpsfwd(trim,fitmodel,vcmmodel,invenu,gps,outdir)
%========================================================
% function[gpsfit]=gpsfwd(trim,fitmodel,vcmmodel,invenu,gps,outdir)
%
%  forward calculation for gps velocities
%
% INPUT:
%  trim:     triangular mesh
%  fitmodel: fitted velocity field
%  vcmmodel: vcm of fitted velocity field
%  invenu:   inversion parameters
%  gps:      gps data
%  outdir:   output directory (optional)
% 
% OUTPUT:
%  gpsfit: fitted gps velocity
%
% Hua Wang @ Uni Leeds, 09/11/2009
%
% 17/03/2015 HW: support multiple gps files
%========================================================
nvtx=length(trim.x);
invs=sum(invenu);
vel = reshape(fitmodel(1:invs*nvtx),nvtx,invs);
vcm = vcmmodel(1:invs*nvtx,1:invs*nvtx);

%forward calculation for GPS
vel_v=reshape(vel,invs*nvtx,1);
ngf=length(gps);

for igf=1:ngf

  ngps=length(gps(igf).site);
  gpsmat=designgps(trim,gps(igf).site,invenu);
  %gpsmat=gpsmat(:,1:invs*nvtx);
  velfit=full(gpsmat*vel_v);
  velfit=reshape(velfit,ngps,invs);
  vcmfit=full(gpsmat*vcm*gpsmat');

  %output GPS
  igpsfit=gps(igf).site;
  obsdim=sum(gps(igf).invenu);
  fidfit=fopen(char(strcat(outdir,'gpsfit',num2str(igf),'.dat')),'w');
  fidobs=fopen(char(strcat(outdir,'gpsobs',num2str(igf),'.dat')),'w');
  for i=1:ngps
    igpsfit(i).vel=velfit(i,1:obsdim);
    %fitted gps data
    std_ve=sqrt(vcmfit(i,i));
    std_vn=sqrt(vcmfit(i+ngps,i+ngps));
    cov_en=vcmfit(i,i+ngps)/std_ve/std_vn;
    %if invs==2
    if obsdim==2
      sel=[i;i+ngps];
      fprintf(fidfit,'%10.4f %10.4f %10.3f %10.3f %10.3f %10.3f %10.5f %s\n',...
             igpsfit(i).lon,igpsfit(i).lat,velfit(i,1),velfit(i,2),std_ve,std_vn,cov_en,char(igpsfit(i).staid));
    else
      sel=[i;i+ngps;i+2*ngps];
      std_vu=sqrt(vcmfit(i+2*ngps,i+2*ngps));
      cov_eu=vcmfit(i,i+2*ngps)/std_ve/std_vu;
      cov_nu=vcmfit(i+ngps,i+2*ngps)/std_vn/std_vu;
      fprintf(fidfit,'%10.4f %10.4f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.5f %10.5f %10.5f %s\n',...
             igpsfit(i).lon,igpsfit(i).lat,velfit(i,1),velfit(i,2),velfit(i,3),std_ve,std_vn,std_vu,cov_en,cov_eu,cov_nu,char(igpsfit(i).staid));
    end
    igpsfit(i).vcm=vcmfit(sel,sel');
  
    %observed gps data
    std_ve=sqrt(gps(igf).site(i).vcm(1,1));
    std_vn=sqrt(gps(igf).site(i).vcm(2,2));
    cov_en=gps(igf).site(i).vcm(1,2)/std_ve/std_vn;
    if obsdim==2
      fprintf(fidobs,'%10.4f %10.4f %10.3f %10.3f %10.3f %10.3f %10.5f %s\n',...
             igpsfit(i).lon,igpsfit(i).lat,gps(igf).site(i).vel(1),gps(igf).site(i).vel(2),std_ve,std_vn,cov_en,char(igpsfit(i).staid));
    else
      std_vu=sqrt(gps(igf).site(i).vcm(3,3));
      cov_eu=gps(igf).site(i).vcm(1,3)/std_ve/std_vu;
      cov_nu=gps(igf).site(i).vcm(2,3)/std_vn/std_vu;
      fprintf(fidobs,'%10.4f %10.4f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.5f %10.5f %10.5f %s\n',...
             igpsfit(i).lon,igpsfit(i).lat,gps(igf).site(i).vel(1),gps(igf).site(i).vel(2),gps(igf).site(i).vel(3),std_ve,std_vn,std_vu,cov_en,cov_eu,cov_nu,char(igpsfit(i).staid));
    end
  end
  %update gpsfit
  gpsfit(igf).invenu=gps(igf).invenu;
  gpsfit(igf).site=igpsfit;
  fclose(fidfit);
  fclose(fidobs);
end
