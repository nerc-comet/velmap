function [interpmat]=designinterp(trim,insar,talk,invenu)
%================================================
%function [interpmat]=designinterp(trim,insar,talk,invenu)
% 
% Design interpolation matrix for insar data
% 
% Input:
%  trim: triangular mesh structure (tri,x,y)
%  insar: insar data structure
%  talk: display message (default 0)
%  invenu: inverse parameters (default [1 1 1])
%
% Output:
%  interpmat: interpolation design matrix
%
% Hua Wang @ Leeds, 17/09/2009
%
% 01/06/2011 HW: replace stackmap by los for nans
% 27/05/2011 HW: set e/n/u contribution one by one
%================================================
if nargin<3
    talk=0;
end

% if no inversion compoents are given, default to all three (enu)
if nargin<4
    invenu=[1 1 1];
end
invs=sum(invenu);
icomp=cumsum([0 invenu(1:2)]);

ntri=length(trim.tri);    %triangular number
nvtx=length(trim.x);      %vertex number
ninsar=length(insar);     %insar rate map number

%design unit vector matrix
disp('calculating unit vectors');
uvec = [];
for i=1:ninsar
    % calculate enu from los and azi
    uvi = unitvec(insar(i).los,insar(i).azi);
    [rows,cols] = size(insar(i).stackmap);
    n = rows*cols;
    
    % reshape dem to column and use to mask unit vectors
%     demv = reshape(insar(i).dem',n,1);
%     uvi(isnan(demv),:) = NaN;
    
    % append
    uvec=[uvec;uvi];  %e/n/u
end

disp('making interpolation kernel');
totpix=length(uvec);
%kernel=sparse(totpix,nvtx);
interpmat=sparse(totpix,invs*nvtx);

% for each triangle in the mesh ...
for itri = 1:ntri
    
    % report progress
    if (mod(itri,100))==0 && (talk==1)
        fprintf('%04d/%04d triangles\n',itri,ntri);
    end
    
    %vertex coordinates for each triangular
    x=trim.x(trim.tri(itri,:)); %x coordinates of the three vertics
    y=trim.y(trim.tri(itri,:)); %y coordinates of the three vertics
    
    % for each insar data set ...
    for is=1:ninsar
        
        % calculate number of grid points between insar map edges and
        % triangle
        ifghdr=insar(is).ifghdr;
        minx = floor((min(x)-ifghdr.xfirst)/ifghdr.xstep)+1;
        maxx = ceil((max(x)-ifghdr.xfirst)/ifghdr.xstep)+1;
        miny = floor((max(y)-ifghdr.yfirst)/ifghdr.ystep)+1;
        maxy = ceil((min(y)-ifghdr.yfirst)/ifghdr.ystep)+1;
        
        %get a patch to speed up searching
        minxp = max(minx,1);
        maxxp = min(maxx,ifghdr.width);
        minyp = max(miny,1);
        maxyp = min(maxy,ifghdr.length);
        
        %first pixel no. for the is-th stackmap
        irow0=0;
        for js=1:is-1
            irow0=irow0+numel(insar(js).los);
        end
        
        %search when triangular is at least partly within the insar track
        if (minxp<maxxp && minyp<maxyp)
            
            %for each point within the patch, check whether it is within the triangular
            for ixp=minxp:maxxp
                for iyp=minyp:maxyp
                    
                    %          if and(~isnan(insar(is).los(iyp,ixp)),~isnan(insar(is).dem(iyp,ixp)))
                    if ~isnan(insar(is).los(iyp,ixp))
                        pt(1)=ifghdr.xfirst+(ixp-1)*ifghdr.xstep; %longitude
                        pt(2)=ifghdr.yfirst+(iyp-1)*ifghdr.ystep; %latitude
                        tri=[x,y];                                %triangular
                        
                        pos=intri(pt,tri);
                        
                        if pos~=0
                            %interpolation function (England & Molnar, 2005, P5, Fm5-8)
                            N=interpk(pt,tri);
                            %row/col number in the design matrix
                            irow=ixp+(iyp-1)*ifghdr.width+irow0;
                            icol=trim.tri(itri,:);
                            %kernel(irow,icol)=N;
                            
                            % for each component ...
                            for ic=1:3
                                if insar(is).proc.invenu(ic)==1
                                    interpmat(irow,icol+icomp(ic)*nvtx) = N.*uvec(irow,ic);
                                end
                            end
                            
                            %set the pixel as nan in the original rate map
                            %so that this pixel will not be evaluated for another triangular
                            insar(is).los(iyp,ixp)=nan;
                        end
                    end
                end
            end
        end
    end
end

for i=1:ninsar
    %  if nnz(and(~isnan(insar(is).los),~isnan(insar(is).dem)))>0
    if nnz(~isnan(insar(is).los))>0
        error('not all points have been evaluated for triangular search, check insar file i=%d',i)
    end
end

disp('making final interpolation matrix');
%disp(size(interpmat))
%disp(size(uvec))
interpmat(isnan(uvec(:,1)),:)=[];      
%disp(size(interpmat))
%remove nans
