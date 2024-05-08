function [p]=intri(tri)
%============================================================
%function [p]=intri(tri)
% calculate the coordinate of the incenter of a triangular
% 
% Input:
%  tri: coordinates of the three vertex [(x1 y1);(x2 y2);(x3 y3)]
%
% Output:
%  p:  coordinate of the incenter (x,y)
%
%The following algorithm is implemented
%  p(x,y) = [a*(xa,ya)+b*(xb,yb)+c*(xc,yc)]/P;
%  a,b,c: lengths of the opposite sides 
%  P=a+b+c;
%  xa,xb,xc: x coordinate of vertics a-c
%  ya,yb,yc: y coordinate of vertics a-c
%
% Hua Wang @ Leeds, 17/11/2009
%============================================================

a(1) = sqrt((tri(2,1)-tri(3,1))^2+(tri(2,2)-tri(3,2))^2); %a_23
a(2) = sqrt((tri(1,1)-tri(3,1))^2+(tri(1,2)-tri(3,2))^2); %b_13
a(3) = sqrt((tri(1,1)-tri(2,1))^2+(tri(1,2)-tri(2,2))^2); %c_12
P=sum(a);
p=a*tri/P;
