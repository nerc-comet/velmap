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

if nargin<4
  invenu=[1 1 1];
end
invs=sum(invenu);
icomp=cumsum([0 invenu(1:2)]);

ntri=length(trim.tri);    %triangular number
nvtx=length(trim.x);      %vertex number
ninsar=length(insar);     %insar rate map number

%design unit verctor matrix
disp('calculating unit vectors');
uvec=[];
for i=1:ninsar
%     if insar(i).proc.atmdegree==1
%     insar(i).los(isnan(insar(i).dem)) = NaN;
%     end
%    disp(insar(i))
%    disp(insar(i))
  [uvi]=unitvec(insar(i).los,insar(i).azi);
  uvec=[uvec;uvi];  %e/n/u
end
%disp(size(uvec))

disp('making interpolation kernel');
totpix=length(uvec);

% create cell arrays for per process insar and interpmat data:
insars = cell(ninsar, 1);
interpmats = cell(ninsar, 1);

% alllocation size for interpmat sparse matrix.
% 0.0005 * size of matrix:
spalloc_size = round(0.0005 * totpix * (invs*nvtx));

% create per process insar and interpmat data:
for i = 1:ninsar
  interpmats{i} = spalloc(totpix, invs*nvtx, spalloc_size);
  insars{i} = insar;
end

% progress message:
disp(['progress (', num2str(ninsar), ' images):']);
fprintf([repmat('.', 1, ninsar), '\n\n']);

% loop ninsar in parallel:
parfor is=1:ninsar

%disp(size(insars{is}(is).los))
%disp(size(insars{is}(is).dem))

  % init irow0 value:
  irow0 = nan;
  % get ifgheader:
  ifghdr=insars{is}(is).ifghdr;
  % get ifg values:
  xfirst = ifghdr.xfirst;
  xstep = ifghdr.xstep;
  ystep = ifghdr.ystep;
  yfirst = ifghdr.yfirst;
  ifg_width = ifghdr.width;

  % loop ntri:
  for itri = 1:ntri

    if (mod(itri,100))==0 && (talk==1)
      fprintf('%04d/%04d triangles\n',itri,ntri);
    end

    % get trim values for this loop:
    trim_itri = trim.tri(itri, :);

    % vertex coordinates for each triangular
    % x coordinates of the three vertics:
    x=trim.x(trim_itri);
    minx=floor((min(x)-xfirst)/xstep)+1;
    maxx=ceil((max(x)-xfirst)/xstep)+1;
    %get a patch to speed up searching
    minxp = max(minx,1);
    maxxp = min(maxx,ifg_width);

    % search when triangular is at least partly within the insar track
    % check x first:
    if (minxp<maxxp)

      % y coordinates of the three vertics
      y=trim.y(trim_itri);
      miny=floor((max(y)-yfirst)/ystep)+1;
      maxy=ceil((min(y)-yfirst)/ystep)+1;
      minyp = max(miny,1);
      maxyp = min(maxy,ifghdr.length);

      % check y:
      if (minyp<maxyp)

        %for each point within the patch, check whether it is within the triangular
        for ixp=minxp:maxxp
          for iyp=minyp:maxyp

            % get iyp -1 value:
            iyp_minus = iyp - 1;

            if ~isnan(insars{is}(is).los(iyp,ixp))
              pt = [ ...
                     xfirst+(ixp-1)*xstep,    ... %longitude
                     yfirst+(iyp_minus)*ystep ... %latitude
                   ];
              tri=[x,y];                          %triangular

              pos=intri(pt,tri);

              if pos~=0

                % get row0 value, if we don't have it:
                if (isnan(irow0))
                  irow0=0;
                  for js=1:is-1
                    irow0=irow0+numel(insars{is}(js).los);
                  end
                end

                %interpolation function (England & Molnar, 2005, P5, Fm5-8)
                N=interpk(pt,tri);
                %row/col number in the design matrix
                irow=ixp+(iyp_minus)*ifg_width+irow0;
                icol=trim_itri;

                %east component
                for ic=1:3
                  if insars{is}(is).proc.invenu(ic)==1
                    interpmats{is}(irow,icol+icomp(ic)*nvtx) = N.*uvec(irow,ic);
                  end
                end

                %set the pixel as nan in the original rate map
                %so that this pixel will not be evaluated for another triangular
                insars{is}(is).los(iyp,ixp)=nan;

              % end 'if pos~=0':
              end
            % end 'if ~isnan(insar(is).los(iyp,ixp))':
            end
          % end 'for iyp=minyp:maxyp':
          end
        % end 'for ixp=minxp:maxxp':
        end
      % end 'if (minyp<maxyp)':
      end
    % end 'if (minxp<maxxp)':
    end
  % end 'for itri = 1:ntri':
  end
  % update progress:
  fprintf('\b|\n');
% end 'parfor is=1:ninsar':
end

% update insar data with data from each process:
for i = 1:ninsar
  insar(i) = insars{i}(i);
end
clear insars;

% update interpmat with data from each process:
interpmat = interpmats{1};
for i = 2:ninsar
  interpmat = interpmat + interpmats{i};
end
clear interpmats;

for i=1:ninsar
  if nnz(~isnan(insar(i).los))>0
    error('not all points have been evaluated for triangular search, check insar file i=%d',i)
  end
end

disp('making final interpolation matrix');
% disp(size(interpmat))
% disp(size(uvec))
% disp(sum(isnan(uvec(:,1))))
% disp(sum(~isnan(uvec(:,1))))
interpmat(isnan(uvec(:,1)),:)=[];            %remove nans
%disp(size(interpmat))