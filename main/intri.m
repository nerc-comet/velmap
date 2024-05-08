function [pos]=intri(p,tri)
%============================================================
%function [pos]=intri(p,tri)
% Check the position of point p w.r.t triangle p1-3
% 
% Input:
%  p:  point (x,y)
%  tri: coordinates of the three vertex [(x1 y1);(x2 y2);(x3 y3)]
%
% Output:
%  pos: position of point P
%     0: outside
%     1: vertex 1
%     2: vertex 2
%     3: vertex 3
%     12: on edge p12
%     23: on edge p23
%     31: on edge p31
%     123: inside
%
%Example:
%pos=intri([0.5 0.5],[0 0; 0 2; 2 0]);
%
%The following algorithm is implemented
% If P is ON or INSIDE the triangle
%      Area(PP1P2) + Area(PP2P3) + Area(PP3P1) = Area(P1P2P3)
% If P is OUTSIDE
%      Area(PP1P2) + Area(PP2P3) + Area(PP3P1) > Area(P1P2P3)
% 
% Area of a triangle can be found using the determinant:
% 
%                        |x1  y1  1|  
%        S = abs ( 1/2 * |x2  y2  1| )
%                        |x3  y3  1|  
%
% Hua Wang @ Leeds, 17/09/2009
%============================================================

s0  = 1/2. *abs(det([tri ones(3,1)]));
if s0==0
  error('input wrong triangle');
end

s12 = 1/2. *abs(det([p  1; tri(1,:) 1; tri(2,:) 1]));
s23 = 1/2. *abs(det([p  1; tri(2,:) 1; tri(3,:) 1]));
s31 = 1/2. *abs(det([p  1; tri(3,:) 1; tri(1,:) 1]));
s=s12+s23+s31;

if s-s0<eps('single')
  pos=1;               %inside 
else
  pos=0;
end

%if s-s0>eps('single')
%  pos=0;               %outside 
%else
%  if s12+s31==0        %vertex 1
%    pos=1;
%  elseif s12+s23==0    %vertex 2
%    pos=2;
%  elseif s23+s31==0    %vertex 3
%    pos=3;
%  elseif s12==0        %on edge p12
%    pos=12;
%  elseif s23==0        %on edge p23
%    pos=23;
%  elseif s31==0        %on edge p31
%    pos=31;
%  else                 %inside
%    pos=123;
%  end
%end
