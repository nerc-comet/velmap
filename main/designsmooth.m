function[smmat] = designsmooth(trim,boundsm,invenu)
%====================================================================
%function[smmat] = designsmooth(trim,boundsm,invenu)
%                                                                    
% Design smoothing matrix for a triangular mesh using Laplacian smoothing 
%
%INPUT
% trim: triangular mesh
% boundsm: boundary smoothing (optional, 0: not smoothing; 1: smoothing (default))
% invenu: inverse parameters (default: [1 1 1])
%OUTPUT
% smmat: smoothing matrix
% 
%see Desbrun et al., SIGGRAPH, 1999, page 321 formula 11 for details
%
%Hua WANG, 28/09/2009
%
% 27/03/2010 HW: make boundary smoothing as optional
%====================================================================

if nargin<2
  boundsm=1;
end
if nargin<3
  invenu=[1 1 1];
end

ntri=size(trim.tri,1);
nvtx=length(trim.x);

%initial smooth matrix
sm = zeros(nvtx,nvtx);

%boundary identify
if boundsm==0
  bound = zeros(nvtx,nvtx,'int8');
end

%make the distance matrix (symmetric matrix)
%note: most of the edges have been used twice, but it will not affect the matrix
index=[1,2,3,1];
for i=1:ntri
  for j=1:3
    %calculate edge distance
    vtx1=trim.tri(i,j);
    vtx2=trim.tri(i,index(j+1));
    sm(vtx1,vtx2)=sqrt((trim.x(vtx1)-trim.x(vtx2))^2+(trim.y(vtx1)-trim.y(vtx2))^2);
    %identify boundary by the calculation times of its distance
    if boundsm==0
      bound(vtx1,vtx2)=bound(vtx1,vtx2)+1;
      bound(vtx2,vtx1)=bound(vtx2,vtx1)+1;
    end
  end
end

%make smoothing factors for each nodes
E=repmat(sum(sm,2),1,nvtx); %sum of the edges for each vertex, E in Desbrun et al., 1999
nz=find(sm);                %non-zeros, diagonal components are still zeros here
sm(nz)=2./sm(nz)./E(nz);    %off-diagonal components, (2/E)*(1/e_ij)
smdiag=-sum(sm,2);          %diagonal components, -2/E*sum(1/e_ij) = -sum(2/E/e_ij)=-sum(sm,2);
sm=sm+diag(smdiag);         %combine off-diagonal and diagonal components

%remove boundaries
if boundsm==0
  %bound==1: the vertex is on the boundary
  boundvtx=sum(bound==1,2);
  sm(boundvtx>0,:)=[];
end

%make the final smoothing matrix for E/N/U components
ncomp=sum(invenu);
smmat=sm;
for i=2:ncomp
  smmat=blkdiag(smmat,sm);
end
%smmat=blkdiag(sm,sm,sm);
