function [N,a,b,c]=interpk(p,tri)
%============================================================
%function [kernel]=interpk(p,tri)
% interpolation kernel for a point p w.r.t triangle p1-3
% 
% Input:
%  p:  point (x,y)
%  tri: coordinates of the three vertex [(x1 y1);(x2 y2);(x3 y3)]
%
% Output:
%  N: interpolation kernel
%  a/b/c: interpolation coefficients
%
% see England & Molnar, 2005, JGR, Page 5, Fms 5-8 for details
% Hua Wang @ Leeds, 28/10/2009
%============================================================

x=tri(:,1);
y=tri(:,2);

a(1)=x(2)*y(3)-x(3)*y(2);
a(2)=x(3)*y(1)-x(1)*y(3);
a(3)=x(1)*y(2)-x(2)*y(1);
b(1)=y(2)-y(3);
b(2)=y(3)-y(1);
b(3)=y(1)-y(2);
c(1)=x(3)-x(2);
c(2)=x(1)-x(3);
c(3)=x(2)-x(1);
delta=x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2));
a=a/delta;
b=b/delta;
c=c/delta;
%N = a+b*lon+c*lat;
N = a+b*p(1)+c*p(2);
